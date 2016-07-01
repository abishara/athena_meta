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

    def __init__(self, options_path, **kwdargs):
        self.options_path = options_path
        self._output_dir = os.path.dirname(self.options_path)

        # set required attributes to None
        for opt in self.required:
          setattr(self, opt, None)

        # set optional to default
        for opt, val in self.optional:
          setattr(self, opt, val)

        self.cluster_settings = ClusterSettings()

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
    def output_dir(self):
        return self._output_dir
    
    @property
    def results_dir(self):
        return os.path.join(self.output_dir, "results")

    @property
    def working_dir(self):
        return os.path.join(self.output_dir, "working")

    @property
    def log_dir(self):
        return os.path.join(self.output_dir, "logs")

    @abc.abstractmethod
    def __getstate__(self):
        """
        allows pickling of Options instances, necessary for ipyparallel
        """
        state = self.__dict__.copy()
        return state

#--------------------------------------------------------------------------
# ref assembly options
#--------------------------------------------------------------------------
class RefAsmOptions(Options):
      
    @classproperty
    def pipe_type(self): return 'ref-asm'

    @classproperty
    def required(self):
      return [
        'longranger_bam_path',
        'longranger_vcf_path',
        'longranger_fqs_path',
        'ref_fasta',
      ]

    @classproperty
    def optional(self):
      return [
        ('regions_bed_path', None),
        ('genome_step_size',   50000),
        ('genome_window_size', 100000),
      ]

    def __init__(self, options_path, debug=False):

        super(RefAsmOptions, self).__init__(options_path)

        self.binaries = None
        self._regions = None

    @property
    def regions(self):
        if self._regions == None and self.regions_bed_path:
          self._regions = util.load_bed(self.regions_bed_path)
        return self._regions

    @property
    def bins_pickle_path(self): 
        return os.path.join(self.working_dir, 'bins.p')

    @property
    def groups_pickle_path(self): 
        return os.path.join(self.working_dir, 'groups.p')
    @property
    def groups2_pickle_path(self): 
        return os.path.join(self.working_dir, 'groups2.p')
    
    def get_group_dir(self, gid, final=False):
      return os.path.join(
        self.working_dir if not final else self.results_dir,
        'groups',
        'group.{}'.format(gid),
      )

    def get_group_fq_dir(self, gid):
      return os.path.join(self.get_group_dir(gid), 'fqs')
    def get_group_asm_dir(self, gid):
      return os.path.join(self.get_group_dir(gid), 'asm')

    def get_bin_dir(self, binid, final=False):
      ctg, b, e, cidx = binid
      return os.path.join(
        self.working_dir if not final else self.results_dir,
        'bins',
        '{}.{}-{}.c{}'.format(ctg, b, e, cidx),
      )
      
    def get_bin_fq_dir(self, binid):
      return os.path.join(self.get_bin_dir(binid), 'fqs')
    def get_bin_asm_dir(self, binid):
      return os.path.join(self.get_bin_dir(binid), 'asm')

    def binary(self, name):
        """
        Checks to see if a path has been specified for an external binary,
        otherwise just return the name of the binary to try running it
        if it's in $PATH
        """
        if self.binaries is not None:
            return self.binaries.get(name, name)
        return name

    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, mode=True):
        self._reference = None
        self._debug = mode

    def __str__(self):
        d = self.serialize()
        d["debug"] = self.debug
        return json.dumps(d, sort_keys=True, indent=4)

    def __getstate__(self):
        """
        allows pickling of Options instances, necessary for ipyparallel
        """
        state = self.__dict__.copy()
        state["_regions"] = None

        return state

#--------------------------------------------------------------------------
# reads options
#--------------------------------------------------------------------------
class ReadsOptions(Options):

    @classproperty
    def pipe_type(self): return 'reads'

    @classproperty
    def required(self):
      return ['tenxfq_path']

    @classproperty
    def optional(self):
      return []

    #def __init__(self, options_path, debug=False):
    #    super(RefAsmOptions, self).__init__(options_path, debug)

    def get_bin_dir(self, binid, final=False):
      assert binid == None
      return os.path.join(
        self.working_dir if not final else self.results_dir,
        'bin',
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
        'bcode_ctg_hits_path',
        'reads_ctg_bam_path',
        'longranger_fqs_path',
      ]

    @classproperty
    def optional(self):
      return []

    #def __init__(self, options_path, debug=False):
    #    super(RefAsmOptions, self).__init__(options_path, debug)

    @property
    def bins_pickle_path(self): 
        return os.path.join(self.working_dir, 'bins.p')
    @property
    def bins2_pickle_path(self): 
        return os.path.join(self.working_dir, 'bins2.p')

    # FIXME hacks to work with collect reads for now
    @property
    def groups_pickle_path(self): 
        return os.path.join(self.working_dir, 'bins.p')
    def get_group_fq_dir(self, binid):
      return os.path.join(self.get_bin_dir(binid), 'fqs')
    # end hacks

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

