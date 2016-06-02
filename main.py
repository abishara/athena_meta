import os
import json
import sys
import collections
import logging

from athena.mlib import util

from athena import pipeline
from athena.options import Options, ClusterSettings

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

    return stages

def clean_up():
  junk = filter(
    lambda(f): (
      f.startswith('sge_controller') or 
      f.startswith('sge_engine') or 
      f.startswith('bcbio-e') or
      f.startswith('bcbio-c')
    ),
    os.listdir('.'),
  )
  map(lambda(f): os.remove(f), junk)

def main(argv):
  """
  1. process command-line arguments
  3. run
  """
  assert len(argv) > 1
  scratch_path = argv[1]
  options = Options(scratch_path, None)
  clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main(sys.argv)

