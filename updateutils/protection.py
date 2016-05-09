"""
This class offers protection against interrupting the update process.

It will generate locks when doing certain processes that will execute a
'cleanup' function before releasing the locks and exiting. Prevents a tree from
becoming polluted to the best of my ability
"""
import sys
from updateutils.exceptions import LockWarning
import warnings

class _Lock(object):
   """ A safety lock """
   def cleanup(self):
      """
      Holder method to avoid any NameError's. This should be the function
      necessary to run in order to clean up this lock
      """
      pass

class GlobalState(object):
   """
   The global state. Protection for certain critical actions should be triggered
   by registering that action with the main GlobalState instance to apply a
   lock.

   Each lock will have associated with it a particular clean-up action to
   perform before the lock can be released.
   
   The GlobalState can only have one lock active at a time. You must release the
   lock (via the release_lock method) when appropriate

   Each lock you want to register needs a method attached to the GlobalState
   class to 'start' that lock. It must be passed the 'cleanup' function that the
   ensuing _Lock must perform
   """
   def __init__(self):
      # We start with no locks
      self.lock = None

   def release_lock(self):
      """ Release the existing lock """
      self.lock = None

   def _lock_decorator(fcn):
      """ Make sure we don't add multiple locks """
      def check_locks(self, *args):
         if self.lock is not None:
            warnings.warn('Lock already engaged. Ignoring the new one.',
                          LockWarning)
            return
         fcn(self, *args)
      return check_locks

   @_lock_decorator
   def test_lock(self, method):
      """ Tests the locking mechanism """
      self.lock = _Lock()
      self.lock.cleanup = method


   @_lock_decorator
   def add_lock(self, method):
      """ Default lock adding mechanism with a cleanup function passed """
      self.lock = _Lock()
      self.lock.cleanup = method

   def cleanup(self):
      """ Cleans up after an interruption or something """
      if self.lock is not None:
         self.lock.cleanup()
         self.lock = None

# Define the global state. This should be used during the main execution to
# register safety locks
global_state = GlobalState()


def signal_handler(sig, frame):
   """ Signal handler to safeguard program termination """
   from threading import Thread
   global global_state

   sys.stderr.write('Caught interrupt signal. Cleaning up\n')

   if global_state.lock is not None:
      # Clean up via a 'thread' to prevent interruptions
      thr = Thread(target=global_state.cleanup)
      thr.start()
      thr.join()

   sys.stderr.write('Exiting\n')
   sys.exit(1)

def test():
   """ Test suite """
   from signal import signal, SIGINT
   from time import sleep
   signal(SIGINT, signal_handler)

   def test_lock():
      sys.stderr.write('Cleaning up test lock. This takes exactly 3 seconds\n')
      sleep(3)
      return

   print 'Instituting the lock in 5 seconds'
   sleep(5)
   global_state.test_lock(test_lock)
   print 'Locked'
   print 'Sleeping for 10 seconds'
   sleep(10)
   global_state.release_lock()
   print 'Lock released!'
   print 'Done with test.'
