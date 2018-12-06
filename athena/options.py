import abc
import os
import json

from athena.mlib import util

class classproperty(object):
  def __init__(self, f):
    self.f = f
  def __get__(self, obj, owner):
    return self.f(owner)

class ClusterSettings(object):
  def __init__(self):
    # local
    self.cluster_type = "local"
    self.processes = 1
    self.cluster_options = {}
  
  @staticmethod
  def local():
    return ClusterSettings()

  @staticmethod
  def multiprocessing(processes):
    settings = ClusterSettings()
    settings.cluster_type = 'multiprocessing'
    settings.processes = processes
    settings.cluster_options = {}
    return settings
    
  @staticmethod
  def deserialize(options_dict):
    settings = ClusterSettings()
    
    if "processes" in options_dict:
      settings.processes = options_dict["processes"]
    if "cluster_type" in options_dict:
      settings.cluster_type = options_dict["cluster_type"]
    if "cluster_options" in options_dict:
      settings.cluster_options = options_dict["cluster_options"]
    
    return settings

#--------------------------------------------------------------------------
# options base
#--------------------------------------------------------------------------
class Options(object):
  __metaclass__ = abc.ABCMeta
  
  @classproperty
  def pipe_type(self): return None
  
  @classproperty
  def required(self):
    """ options required to be specified """
    raise Exception("Not implemented yet")
  
  @classproperty
  def optional(self):
    """ options that are optional to specify. tuple of (key, default
    value) """
    raise Exception("Not implemented yet")

  def __init__(self, options_path,  **kwdargs):

    self.options_path = options_path
    self._output_dir = os.path.dirname(self.options_path)
    if self._output_dir == '':
      self._output_dir = './'
    
    # set required attributes to None
    for opt in self.required:
      setattr(self, opt, None)
    
    # set optional to default
    for opt, val in self.optional:
      setattr(self, opt, val)

    # set any arguments specified by keyword args
    for opt, val in kwdargs.items():
      assert opt in self.required or opt in self.optional, \
        "unexpected {} keyword argument {}".format(
          type(self).__name__, opt,
        )
      setattr(self, opt, val)
    
    self.cluster_settings = ClusterSettings()

  def set_cl_args(self, args):
    self.force_reads = args.force_reads
    if args.threads > 1:
      self.cluster_settings = ClusterSettings.multiprocessing(args.threads)

  @classmethod
  def deserialize(cls, options_path):
    # load json config
    with open(options_path) as f:
      options_dict = json.load(f)
  
    options = cls(options_path)
    # required
    for opt in cls.required:
      assert opt in options_dict, 'required option "{}" missing'.format(opt)
      setattr(options, opt, options_dict[opt])
  
    # optional
    for opt, val in cls.optional:
      setattr(options, opt, options_dict.get(opt, val))
  
    # cluster settings
    options.cluster_settings = ClusterSettings.deserialize(
      options_dict.get("cluster_settings", {}))
  
    return options

  @property
  def output_dir(self): return self._output_dir
  
  @property
  def results_dir(self): return os.path.join(self.output_dir, "results")
  
  @property
  def working_dir(self): return os.path.join(self.output_dir, "working")
  
  @property
  def log_dir(self): return os.path.join(self.output_dir, "logs")
  
  @abc.abstractmethod
  def __getstate__(self):
    """
    allows pickling of Options instances, necessary for ipyparallel
    """
    state = self.__dict__.copy()
    return state

#--------------------------------------------------------------------------
# metagenome assembly options
#--------------------------------------------------------------------------
class MetaAsmOptions(Options):

  @classproperty
  def pipe_type(self): return 'meta-asm'
  
  @classproperty
  def required(self):
    return [
      'ctgfasta_path',
      'reads_ctg_bam_path',
      'input_fqs',
    ]
  
  @classproperty
  def optional(self):
    return [
      ('cheat_seeds', None),
      ('force_reads', False),
      ('ds_subasm_cov', 100),
      ('seed_self_asm_size', 10000),
    ]
  
  def __init__(self, options_path, **kwdargs):
      super(MetaAsmOptions, self).__init__(options_path, **kwdargs)
  
  @property
  def bins_pickle_path(self): 
      return os.path.join(self.working_dir, 'bins.p')
  
  @property
  def groups_pickle_path(self): 
      return os.path.join(self.working_dir, 'bins.p')
  
  def get_bin_dir(self, binid, final=False):
    assert binid.startswith('bin')
    return os.path.join(
      self.working_dir if not final else self.results_dir,
      binid,
    )
    
  def get_bin_fq_dir(self, binid):
    return os.path.join(self.get_bin_dir(binid), 'fqs')
  def get_bin_asm_dir(self, binid):
    return os.path.join(self.get_bin_dir(binid), 'asm')
  
  def __str__(self):
    return self.__class__.__name__
  
  def __getstate__(self):
    """
    allows pickling of Options instances, necessary for ipyparallel
    """
    state = self.__dict__.copy()
    return state

