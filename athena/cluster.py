import multiprocessing

class Cluster(object):
  def __init__(self, processes, cluster_settings):
    self.processes = processes
    self.cluster_settings = cluster_settings


class LocalCluster(Cluster):
  def map(self, fn, args):
    for arg in args:
      fn(arg)


class IPCluster(Cluster):
  def __init__(self, processes, cluster_settings):
    super(IPCluster, self).__init__(processes, cluster_settings)

  def map(self, fn, args):
    from cluster_helper.cluster import cluster_view

    cluster_args = {
      "scheduler": 'slurm',
      "queue": 'owners',
      #"scheduler": 'sge',
      #"queue": '',
      "num_jobs": self.processes,
      #"extra_params": {"run_local": True}
    }

    cluster_args.update(self.cluster_settings.cluster_options)
    
    print cluster_args
    
    with cluster_view(**cluster_args) as view:
      async_results = view.map(fn, args, block=False)
      async_results.wait_interactive()
      return async_results.get()

class MultiprocessingCluster(Cluster):
  def map(self, fn, args):
    pool = multiprocessing.Pool(processes=self.processes)
    return pool.map_async(fn, args).get(9999999)

