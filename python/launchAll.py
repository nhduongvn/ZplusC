import os

task = 'haddtree'
#task = 'singleprep'
for line in open('datalist.txt').readlines():
  line = line.replace("\n","")
  if line.find("#") != -1:
    continue
  print line
  os.system("python submitThem.py TestZllHbb13TeV " + task + " --sample "  + line + " --verbose -n 10")
