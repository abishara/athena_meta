import os
import sys
import collections
import logging

from athena.mlib import util

from athena import pipeline
from athena.options import RefAsmOptions, ReadsOptions, MetaAsmOptions

from athena.stages import haplotype_reads
from athena.stages import collect_reads
from athena.stages import assemble_bins
from athena.stages import group_bins
from athena.stages import assemble_groups
from athena.stages import bin_meta_reads
from athena.stages import assemble_meta_bins

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

    if options.pipe_type == 'ref-asm':
      stages["haplotype_reads"] = haplotype_reads.HaplotypeReadsStep
      stages["collect_bin_reads"] = collect_reads.CollectBinReadsStep
      stages["assemble_bins"] = assemble_bins.AssembleBinnedStep
      stages["group_bins"] = group_bins.GroupBinsStep
      stages["collect_group_reads"] = collect_reads.CollectGroupReadsStep
      stages["assemble_groups"] = assemble_groups.AssembleGroupsStep
    elif options.pipe_type == 'meta-asm':
      stages["bin_reads"] = bin_meta_reads.BinMetaReadsStep
      stages["collect_bin_reads"] = collect_reads.CollectGroupReadsStep
      stages["assemble_bins"] = assemble_meta_bins.AssembleMetaBinnedStep
    elif options.pipe_type == 'reads':
      stages["assemble_reads"] = assemble_bins.AssembleSpecReadsStep
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
  usage: athena.py <path/to/config.json> [pipeline]

  pipeline: {ref-asm, reads, meta-asm}, default: ref-asm

  NOTE: dirname(config.json) specifies root output directory
  '''
  if len(argv) != 2 and len(argv) != 3:
    print help_str
    sys.exit(1)

  config_path = argv[1]
  pipe_type = 'ref-asm'
  if len(argv) == 3:
    pipe_type = argv[2]
    if pipe_type not in [
      'ref-asm',
      'reads',
      'meta-asm',
    ]:
      print >> sys.stderr, 'error: incorrect pipeline specified'
      sys.exit(2)

  # load config json
  options_cls = {
    'ref-asm' : RefAsmOptions,
    'reads'   : ReadsOptions,
    'meta-asm': MetaAsmOptions,
  }[pipe_type]
  options = options_cls.deserialize(config_path)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main(sys.argv)

