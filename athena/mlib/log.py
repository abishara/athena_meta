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
  
  def log(self, message):
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message_str = "{} - {}\n".format(date_str, message)
    
    self.log_file.write(message_str)
    self.log_file.flush()
    
    sys.stderr.write(message_str)
    sys.stderr.flush()
  
  def exception(self, exception):
    trace_string = traceback.format_exc()
    self.log("="*10 + " Exception " + "="*10)
    for line in trace_string.split("\n"):
        self.log(line)
    
    self.log(exception)
  
  def error(self, error_message):
    self.log("ERROR: {}".format(error_message))

