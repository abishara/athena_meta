import sys
import datetime
import traceback

class Logger(object):
  """
  This replaces the built-in logging module, which doesn't work reliably 
  across ipyparallel jobs
  """
  
  def __init__(self, log_path):
    self.log_path = log_path
    self.log_file = None
    self.open_log()
  
  def open_log(self):
    if self.log_file is None:
      self.log_file = open(self.log_path, "w")

  def close_log(self):
    if self.log_file != None:
      self.log_file.close()
  
  def _fmt_msg(self, msg, lvl, end='\n'):
    assert lvl in set(['info', 'debug', 'error'])
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return "{} - {} - {}{}".format(date_str, lvl.upper(), msg, end)
    
  def _log_write(self, msg):
    self.log_file.write(msg)
    self.log_file.flush()

  def _std_write(self, msg, std=sys.stdout):
    std.write(msg)
    std.flush()

  def broadcast(self, msg, end='\n'):
    msg = self._fmt_msg(msg, 'info', end)
    self._log_write(msg)
    self._std_write(msg)

  def info(self, msg, end='\n'):
    msg = self._fmt_msg(msg, 'info', end)
    self._log_write(msg)
    #self._std_write(msg)

  def debug(self, msg, end='\n'):
    msg = self._fmt_msg(msg, 'debug', end)
    self._log_write(msg)
  
  def error(self, msg, end='\n'):
    msg = self._fmt_msg(msg, 'error', end)
    self._log_write(msg)
    self._std_write(msg, sys.stderr)
  
  def exception(self, exception):
    trace_string = traceback.format_exc()
    self.error("="*10 + " Exception " + "="*10)
    for line in trace_string.split("\n"):
        self.error(line)
    self.error(exception)

