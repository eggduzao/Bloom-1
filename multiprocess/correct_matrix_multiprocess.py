

# Python
import os
import gc
import sys
import codecs
import traceback
import subprocess
import configparser
import multiprocessing


def change(mydict, a, b):

  mydict[a] = b

def add_change(queue, mydict, a, b):

  queue.append((mydict, a, b))

def run_change(ncpu, queue):
    
  # Execute job queue
  pool = multiprocessing.Pool(ncpu)
  dump_process_output = pool.starmap(change, [arguments for arguments in queue])
  pool.close()
  pool.join()


###################################################################################################
# Execution
###################################################################################################

if __name__ == "__main__":

  manager = multiprocessing.Manager()
  mydict = manager.dict()
  queue = []
  ncpu = 4

  for i in range(0,10):
    add_change(queue, mydict, str(i), i)

  run_change(ncpu, queue)

  print(mydict)

