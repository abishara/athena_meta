import os

class ClusterSettings(object):
    def __init__(self):

        self.cluster_type = "local"
        self.processes = 2
        self.cluster_options = {}

        # ipython cluster
        self.cluster_type = "IPCluster"
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
    def __init__(self, scratch_path, options_path, debug=False):
        self.options_path = options_path
        self._output_dir = scratch_path
        #self._output_dir = os.path.dirname(self.options_path)

        self.genome_step_size = 50000
        self.genome_window_size = 10000
        self.samples = {}
        self.ref_fasta = '/srv/gsfs0/projects/batzoglou/abishara/data/reference/hg38/refdata-hg38/fasta/genome.fa'
        self.bwa_index = None
        self.binaries = None
        self._reference = None
        self._constants = None

        self.cluster_settings = ClusterSettings()

        self.debug = debug


    def serialize(self, ):
        samples = dict((name, self.samples[name].serialize())
                       for name in self.samples)

        d = {"samples": samples,
             "ref_fasta": self.ref_fasta,
             "cluster_settings": self.cluster_settings.serialize(),
             "bwa_index": self.bwa_index,
             "binaries": self.binaries
        }

        return d

    @staticmethod
    def deserialize(options_dict, options_path):
        samples = {}
        for sample_name in options_dict["samples"]:
            sample_info = options_dict["samples"][sample_name]
            cursample = Sample(sample_name)
            for dataset in sample_info:
                cursample.datasets.append(
                    svdatasets.Dataset.deserialize(cursample, dataset))

            samples[sample_name] = cursample

        options = Options(options_path)
        options.samples = samples
        options.ref_fasta = options_dict["ref_fasta"]
        options.cluster_settings = ClusterSettings.deserialize(
            options_dict.get("cluster_settings", {}))

        options.bwa_index = options_dict.get("bwa_index", None)
        options.binaries = options_dict.get("binaries", None)

        return options

    def iter_10xdatasets(self):
        for sample_name, sample in self.samples.items():
            for dataset in sample.datasets:
                if isinstance(dataset, svdatasets.TenXDataset):
                    yield sample, dataset
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
    
    def sample_info(self, sample_name):
        """
        Gets some pre-calculated values for the given sample (eg frag length 
        distributions)
        """
        info_path = os.path.join(
            self.results_dir, "sample_info.{}.pickle".format(sample_name))

        if not os.path.exists(info_path):
            raise Exception("SampleInfoStep needs to be run first!")

        return utilities.pickle.load(open(info_path))

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

        return state

