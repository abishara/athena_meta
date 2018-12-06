import os
import sys
import collections
import logging
import argparse

import subprocess
import gzip
import tempfile
import shutil
import pysam

from athena.mlib import util

from athena import pipeline
from athena.options import MetaAsmOptions

from athena.stages import check_reads
from athena.stages import subassemble_reads
from athena.stages import assemble_olc

logging.basicConfig(format='%(message)s', level=logging.DEBUG)

def run(options):
  """
  1. create output directories
  2. collect args for each stage
  3. check which stages need to run
  4. iterate through stages and submit jobs
  5. validate that we're done running
  """
  
  util.mkdir_p(options.output_dir)
  util.mkdir_p(options.results_dir)
  util.mkdir_p(options.working_dir)
  util.mkdir_p(options.log_dir)
  
  stages = get_stages(options)
  runner = pipeline.Runner(options)
  
  for stage_name, stage in stages.items():
    runner.run_stage(stage, stage_name)

def clean(options):
  stages = get_stages(options)

  for stage_name, stage in stages.items():
    stage.clean_all_steps(options)

def get_stages(options):
  stages = collections.OrderedDict()
  
  if options.pipe_type == 'meta-asm':
    stages["check_reads"] = check_reads.CheckReadsStep
    stages["subassemble_reads"] = subassemble_reads.SubassembleReadsStep
    stages["assemble_olc"] = assemble_olc.AssembleOLCStep
  else:
    raise Exception("Pipeline not implemented yet")
  
  return stages

def clean_up():
  junk = filter(
    lambda(f): (
      f.startswith('SLURM_controller') or 
      f.startswith('SLURM_engine') or 
      f.startswith('sge_controller') or 
      f.startswith('sge_engine') or 
      f.startswith('bcbio-')
    ),
    os.listdir('.'),
  )
  map(lambda(f): os.remove(f), junk)

def test(args):
  logging.info('running tiny test assembly')
  src_testdir_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'test_data',
  )
  # set up test directory and copy over test reads from src
  testd = tempfile.mkdtemp(prefix='athena-test')
  fqgz_path = os.path.join(src_testdir_path, 'reads.fq.gz')
  seedgz_path = os.path.join(src_testdir_path, 'seeds.fq.gz')

  for fname in ['reads.fq.gz', 'seeds.fa.gz']:
    in_path = os.path.join(src_testdir_path, fname)
    out_path = os.path.join(testd, fname.strip('.gz'))
    with gzip.open(in_path, 'r') as fin, \
         open(out_path, 'w') as fout:
      fout.write(fin.read())

  # create bwa index of seeds and align reads
  readsfq_path = os.path.join(testd, 'reads.fq')
  seedfa_path = os.path.join(testd, 'seeds.fa')
  bam_path = os.path.join(testd, 'align-reads.seeds.bam')
  cmd = 'bwa index {}'.format(seedfa_path)
  subprocess.check_call(cmd, shell=True)
  cmd = 'bwa mem -p -C {0} {1} | samtools sort -o {2} - && samtools index {2}'.format(
    seedfa_path, readsfq_path, bam_path,
  )
  subprocess.check_call(cmd, shell=True)

  # format options for test directory and run
  options_path = os.path.join(testd, 'dummy.json')
  options = MetaAsmOptions(
    options_path,
    ctgfasta_path=seedfa_path,
    reads_ctg_bam_path=bam_path,
    input_fqs=readsfq_path,
  )
  options.set_cl_args(args)
  try:
    run(options)
  except Exception, e:
    logging.error('='*30)
    logging.error('test failed to run to completion')
    sys.exit(1)

  # check outputs for correctness
  outfa_path = os.path.join(testd, 'results/olc/athena.asm.fa')
  seedlen = sum(pysam.FastaFile(seedfa_path).lengths)
  asmlen = max(pysam.FastaFile(outfa_path).lengths)
  if seedlen > asmlen:
    logging.error('test ran to completion, but failed to assemble seed contigs')
    logging.error('  output run in {}'.format(testd))
    return False
  else:
    logging.info('--> test completed successfully.\n')
    shutil.rmtree(testd)
    return True

#--------------------------------------------------------------------------
# main
#--------------------------------------------------------------------------
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config',
    default=None,
    help='input JSON config file for run, NOTE: dirname(config.json) specifies root output directory'
  )
  parser.add_argument(
    '--check_prereqs',
    action='store_true',
    help='test if external deps visible in environment',
  )
  parser.add_argument(
    '--test',
    action='store_true',
    help='run tiny assembly test to check setup and prereqs'
  )
  parser.add_argument(
    '--force_reads',
     action='store_true',
     help='proceed with subassembly even if input *bam and *fastq do not pass QC'
  )
  parser.add_argument(
    '--threads',
    type=int,
    default=1,
    help='number of multiprocessing threads',
  )
  if len(sys.argv) == 1:
    print 'error: one of either --config, --check_prereqs, or --test must be specified\n'
    parser.print_help(sys.stderr)
    sys.exit(1)

  args = parser.parse_args()
  if sum([args.test, args.check_prereqs, (args.config != None)]) > 1:
    print 'only one of --config, --check_prereqs, --test can be specified\n'
    parser.print_help(sys.stderr)
    sys.exit(1)

  if args.check_prereqs:
    passed = util.check_prereqs()
    if passed:
      sys.exit(0)
    else:
      print 'failed prereq checks, external deps required'
      sys.exit(1)

  if args.test:
    test(args)
    sys.exit(0)

  assert os.path.isfile(args.config), \
    "config {} is not a path to a file".format(args.config)

  # load config json
  try:
    options = MetaAsmOptions.deserialize(args.config)
  except Exception as e:
    print >> sys.stderr, "{} invalid JSON config file, Exception:".format(args.config)
    print >> sys.stderr, str(e)
    sys.exit(2)

  # set command line overrides
  options.set_cl_args(args)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main()

