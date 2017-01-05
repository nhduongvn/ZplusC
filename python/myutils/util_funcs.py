import ROOT
import os,sys
import math

def deltaPhi(phi1, phi2):
    '''result = phi1 - phi2
    if result > ROOT.TMath.Pi():
      result -= 2*ROOT.TMath.Pi()
    elif result <= -ROOT.TMath.Pi():
      result += 2*ROOT.TMath.Pi()
    return result
    '''    
    result = phi1 - phi2
    while (result > math.pi): result -= 2*math.pi
    while (result <= -math.pi): result += 2*math.pi
    return result

def deltaR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = deltaPhi(phi1, phi2)
    return ROOT.TMath.Sqrt(deta*deta + dphi*dphi);

def sortPt(idxs, pts): #idxs are the index in pts array. For example pts are pts of all jets and idxs is array contain selected jet idex in pts array(for example after pt, eta cut)
  for i in range(0, len(idxs)):
    for j in range(i+1, len(idxs)):
      idx_i = idxs[i]
      idx_j = idxs[j] 
      if pts[idx_j] > pts[idx_i]:
        idxs[i] = idx_j
        idxs[j] = idx_i

def sortPtFatjet(idxs, pts, nJet):
  idxTmps = []
  for i in range(0, nJet):
    idxTmps.append(idxs[i])
  sortPt(idxTmps,pts)
  for i in range(0,nJet):
    idxs[i] = idxTmps[i]


def findSubFolders(path, fileList,eos=True):
  if path.find('failed') != -1: return 1
  os.system('rm tmp.txt')
  if eos: os.system('xrdfs root://cmseos.fnal.gov ls -u ' + path + '/ > tmp.txt')
  else: os.system('ls ' + path + '/ > tmp.txt')
  
  f = open('tmp.txt')
  lines = f.readlines()
  f.close()
  #print lines
  if len(lines) == 0:
    return 0
  
  for line in lines:
    line = line.replace('\n','')
    if line.find('log.tar.gz') != -1: 
      continue
    if not eos: line = path + '/' + line
    if line.find('.root') != -1: 
      os.system('echo \'' + line + '\' >> ' + fileList)
    else:
      if eos:
        line = '/store' + line.split('/store')[1]
      findSubFolders(line,fileList)
     
  return 1




'''import random
a = []
for i in range(1,11): a.append(random.uniform(1,100))

print a

idxs = [3,4,8,9]

print "Before: "
for i in range(0,len(idxs)): print a[idxs[i]]

sortPt(idxs,a)

print "After: "

for i in range(0,len(idxs)): print a[idxs[i]]
'''

