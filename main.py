import os
import sys
import collections
import logging

from athena.mlib import util

from athena import pipeline
from athena.options import MetaAsmOptions, MetaHapOptions

from athena.stages import index_reads
from athena.stages import bin_meta_reads
from athena.stages import assemble_meta_bins
from athena.stages import assemble_olc

from athena.stages import haplotype_reads

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
      stages["bin_reads"] = bin_meta_reads.BinMetaReadsStep
      stages["index_reads"] = index_reads.IndexReadsStep
      #stages["assemble_bins"] = assemble_meta_bins.AssembleMetaBinnedStep
      stages["assemble_olc"] = assemble_olc.AssembleOLCStep
    elif options.pipe_type == 'meta-hap':
      stages["call_variants"] = haplotype_reads.CallVariantsStep
      stages["haplotype_reads"] = haplotype_reads.HaplotypeReadsStep
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

def main(argv):
  """
  1. process command-line arguments
  3. run
  """

  help_str = '''
  usage: athena_meta.py <path/to/config.json> [pipeline]

  pipeline: {meta-asm, meta-hap}, default: meta-asm

  NOTE: dirname(config.json) specifies root output directory
  '''
  if len(argv) != 2 and len(argv) != 3:
    print help_str
    sys.exit(1)

  config_path = argv[1]
  pipe_type = 'meta-asm'
  if len(argv) == 3:
    pipe_type = argv[2]
    if pipe_type not in [
      'meta-asm',
      'meta-hap',
    ]:
      print >> sys.stderr, 'error: incorrect pipeline specified'
      sys.exit(2)

  # load config json
  options_cls = {
    'meta-asm': MetaAsmOptions,
    'meta-hap': MetaHapOptions,
  }[pipe_type]
  options = options_cls.deserialize(config_path)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main(sys.argv)

