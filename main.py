import os
import sys
import collections
import logging

from athena.mlib import util

from athena import pipeline
from athena.options import Options, ClusterSettings

from athena.stages import haplotype_reads
from athena.stages import collect_reads
from athena.stages import assemble_bins

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

    stages = get_stages()
    runner = pipeline.Runner(options)

    for stage_name, stage in stages.items():
        runner.run_stage(stage, stage_name)

def clean(options):
  stages = get_stages()

  for stage_name, stage in stages.items():
    stage.clean_all_steps(options)


def get_stages():
    stages = collections.OrderedDict()

    stages["haplotype_reads"] = haplotype_reads.HaplotypeReadsStep
    stages["collect_reads"] = collect_reads.CollectReadsStep
    stages["assemble_bins"] = assemble_bins.AssembleBinsStep

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
  usage: athena.py <path/to/config.json>

  NOTE: dirname(config.json) specifies root output directory
  '''
  if len(argv) != 2:
    print help_str
    sys.exit(1)

  config_path = argv[1]

  # load config json
  options = Options.deserialize(config_path)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main(sys.argv)

