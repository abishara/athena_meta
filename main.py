import os
import sys
import collections
import logging

import subprocess
import gzip
import tempfile
import shutil
import pysam

from athena.mlib import util

from athena import pipeline
from athena.options import MetaAsmOptions

from athena.stages import index_reads
from athena.stages import bin_meta_reads
from athena.stages import assemble_meta_bins
from athena.stages import assemble_olc

logging.basicConfig(format='%(message)s', level=logging.DEBUG)

def test():
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
    logging.error('test failed to assemble seed contigs')
    logging.error('  output run in {}'.format(testd))
  else:
    logging.info('--> test completed successfully.\n')

  #print 'ran in testd', testd
  shutil.rmtree(testd)

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
    stages["bin_reads"] = bin_meta_reads.BinMetaReadsStep
    stages["index_reads"] = index_reads.IndexReadsStep
    stages["assemble_bins"] = assemble_meta_bins.AssembleMetaBinnedStep
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

def main():
  """
  1. process command-line arguments
  3. run
  """
  #test()

  argv = sys.argv
  help_str = '''
  usage: athena-meta <path/to/config.json>

  NOTE: dirname(config.json) specifies root output directory
  '''
  if len(argv) != 2:
    print help_str
    sys.exit(1)

  config_path = argv[1]
  if not os.path.isfile(config_path):
    print >> sys.stderr, "must specify valid config_path"
    print >> sys.stderr, help_str
    sys.exit(1)

  # load config json
  options_cls = {
    'meta-asm': MetaAsmOptions,
  }['meta-asm']
  try:
    options = options_cls.deserialize(config_path)
  except Exception as e:
    print >> sys.stderr, "{} invalid JSON config file, Exception:".format(config_path)
    print >> sys.stderr, str(e)
    print >> sys.stderr, help_str
    sys.exit(2)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main()

