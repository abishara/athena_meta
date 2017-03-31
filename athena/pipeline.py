import logging

from athena import cluster

def make_cluster(cluster_settings, processes):
  cluster_classes = {
    "IPCluster": cluster.IPCluster,
    "local": cluster.LocalCluster,
    "multiprocessing": cluster.MultiprocessingCluster,
  }
  
  return cluster_classes[cluster_settings.cluster_type](
      processes, cluster_settings)

class Runner(object):
  def __init__(self, options):
    self.options = options
    self.cluster = None
  
  def run_stage(self, stage, stage_name):
    to_run = []
    for step_chunk in stage.get_steps(self.options):
      if step_chunk.needs_to_run():
        to_run.append(step_chunk)
    
    logging.info("{} {} {}".format("="*30, stage_name, "="*30))
    
    if len(to_run) > 0:
      logging.info("{} chunks to run. Starting...".format(len(to_run)))
      
      processes = min(self.options.cluster_settings.processes, len(to_run))
      cluster = make_cluster(self.options.cluster_settings, processes)
      cluster.map(_run_chunk, to_run)
      
      logging.info("--> {} completed.\n".format(stage_name))
    else:
      logging.info("--> 0 chunks need to be run. Skipping...\n")

    if stage.deliver_message(self.options):
      logging.info(stage.deliver_message(self.options))

def _run_chunk(chunk):
  try:
    chunk.start_logging()
  except Exception, e:
    print "== Error starting logging =="
    raise
  
  try:
    chunk.run()
    chunk.finalize()
    chunk.stop_logging()
  except Exception, e:
    chunk.logger.exception(e)
    raise
      
  if chunk.needs_to_run():
    chunk.logger.error("Step {} failed to produce output.".format(chunk))

