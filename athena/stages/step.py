import abc
import datetime
import os
import shutil
import time

from ..mlib import log


class StepChunk(object):
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def get_steps(options):
        """ """
        raise Exception("this abstract staticmethod needs to be instantiated by subclasses")

    @abc.abstractmethod
    def __init__(self, options, **kwdargs):
        """ must take as input the options instance and a dict of arguments """
        pass

    @abc.abstractmethod
    def outpaths(self, final):
        """ return a dictionary of names -> output paths; final indicates if the
        paths should be for the temporary, working version (final=False) or the
        completed final version (final=True) """
        return

    @abc.abstractmethod
    def run(self):
        """ actually runs the pipeline step on the arguments, given the options, 
        defined in __init__ """
        pass

    def __str__(self):
        raise Exception("not implemented")

    @property
    def log_path(self):
        return os.path.join(self.options.log_dir, str(self))

    def start_logging(self):
        self.logger = log.Logger(self.log_path)
        self.logger.log("--starting logging {} --".format(str(self)))
        self._start_time = time.time()

    def stop_logging(self):
        elapsed = time.time() - self._start_time
        self.logger.log("-> finished running step; time elapsed: {}".format(datetime.timedelta(seconds=elapsed)))
        self.logger.log("--stopping logging--")
        

    @classmethod
    def clean_all_steps(cls, options):
        for chunk in cls.get_steps(options):
            chunk.clean()
            
    def clean(self):
        import logging
        for path_name, path in self.outpaths(final=True).items():
            if os.path.exists(path):
                logging.info("Removing {}".format(path))
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)

        if os.path.exists(self.log_path):
            os.remove(self.log_path)

    def needs_to_run(self):
        """ checks if any of the output files are missing """
        paths = self.outpaths(final=True)
        for name, path in paths.items():
            if not os.path.exists(path):
                return True

        return False

    def finalize(self):
        temp_paths = self.outpaths(final=False)
        final_paths = self.outpaths(final=True)

        for key, working_path in temp_paths.items():
            if working_path != final_paths[key]:
                if not os.path.exists(working_path):
                    raise Exception("{} step failed to produce output: {}".format(
                        self.__class__.__name__, working_path))
                shutil.move(working_path, final_paths[key])
                time.sleep(0.1)

        assert not self.needs_to_run()
