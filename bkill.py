import os
import subprocess

list = subprocess.check_output(['bjobs'])
list = list.split('\n')

for line in list:
  ID = line.split(' ')[0]

  if not ID.isdigit():
    continue
  
  cmd = 'bkill '+ID
  print cmd
  os.system(cmd)
