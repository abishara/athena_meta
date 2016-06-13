import os
import json

from athena.mlib import util

class ClusterSettings(object):
    def __init__(self):

        # mp cluster
        self.cluster_type = "multiprocessing"
        self.processes = 2
        self.cluster_options = {}

        # ipython cluster
        self.cluster_type = "IPCluster"
        self.processes = 2
        self.processes = 70
        self.cluster_options = {}

        # local
        self.cluster_type = "local"
        self.processes = 2
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

    def serialize(self):
        return {
            "processes": self.processes,
            "cluster_type": self.cluster_type,
            "cluster_options": self.cluster_options
        }


class Options(object):
      
    def __init__(self, options_path, debug=False):
        self.options_path = options_path
        self._output_dir = os.path.dirname(self.options_path)

        # inputs from longranger
        self.longranger_bam_path = None
        self.longranger_vcf_path = None
        self.longranger_fqs_path = None

        self.regions_bed_path = None

        self.genome_step_size = 50000
        self.genome_window_size = 100000
        self.ref_fasta = None
        self.binaries = None
        self._reference = None
        self._constants = None

        self.cluster_settings = ClusterSettings()

        self.debug = debug

        self._regions = None

    def serialize(self, ):
        d = {
          "ref_fasta": self.ref_fasta,
          "cluster_settings": self.cluster_settings.serialize(),
          "binaries": self.binaries
        }

        return d

    @staticmethod
    def deserialize(options_path):
        # load json config
        with open(options_path) as f:
          options_dict = json.load(f)

        options = Options(options_path)
        # required
        options.ref_fasta = options_dict["ref_fasta"]
        options.longranger_bam_path = options_dict["longranger_bam_path"]
        options.longranger_vcf_path = options_dict["longranger_vcf_path"]
        options.longranger_fqs_path = options_dict["longranger_fqs_path"]
        # optional
        options.regions_bed_path = options_dict.get("regions_bed_path", None)
        options.binaries = options_dict.get("binaries", None)

        options.cluster_settings = ClusterSettings.deserialize(
            options_dict.get("cluster_settings", {}))


        return options

    @property
    def regions(self):
        if self._regions == None and self.regions_bed_path:
          self._regions = util.load_bed(self.regions_bed_path)
        return self._regions

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
    def bins_pickle_path(self): 
        return os.path.join(self.working_dir, 'bins.p')


    @property
    def log_dir(self):
        return os.path.join(self.output_dir, "logs")

    @property
    def reference(self):
        if self._reference is None:
            self._reference = Reference(self.ref_fasta, self.debug)
        return self._reference
    
    @property
    def constants(self):
        if self._constants is None:
            self._constants = constants.get_constants(self)
        return self._constants
    
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
        state["_reference"] = None
        state["_constants"] = None
        state["_regions"] = None

        return state

