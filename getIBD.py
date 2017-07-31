#!/usr/bin/python
import sys
import argparse
import os.path


parser=argparse.ArgumentParser(
    description='''getIBD.py, 7-11-2017''',
    epilog="""Requires the output of three Refined IBD runs (using different seeds on same datasets) and a map file of all markers used in the Refined IBD runs.\nCreates output file *.ibd12 which contains columns sample1, sample2, IBD1 length/proportion, IBD2 length/proportion, and optionally creates a *.seg which contains sample1, sample2, chromosome, IBD type, region start (in cM), region end (in CM).""")
parser.add_argument('-o', type=str, nargs=1, default=['out'], help='Output file prefix', metavar='out')
parser.add_argument('-f', type=str, nargs=3, required=True, help='Three Refined IBD files', metavar=('file1','file2','file3'))
parser.add_argument('-m', type=str, nargs=1, required=True, help='Map file (PLINK format), non-zero cM column required', metavar='file.map')
parser.add_argument('-t', type=int, nargs=1, default=[0], help='Whether to print IBD1 & IBD2 lengths in cM (=0), or IBD1 & IBD2 proportions (=1)',metavar='0/1')
parser.add_argument('-s', type=int, nargs=1, default=[0], help='Whether to print *.seg file which contains the discovered pairwise IBD regions and whether they are IBD1 or IBD2', metavar='0/1')
#args=parser.parse_args(['-f','./run1/test.ibd', './run2/test.ibd', './run3/test.ibd','-m','safs_filter2_geno0.02_mind0.1_ALL_fixed.map'])
args=parser.parse_args()

print('\n')
if os.path.isfile(args.o[0]+'.ibd12') or os.path.isfile(args.o[0]+'.seg'):
  k=1
  new = args.o[0] + '-' + str(k)
  while os.path.isfile(new+'.ibd12') or os.path.isfile(new+'.seg'):
    k = k + 1
    new = args.o[0] + '-' + str(k)
  print('Output file with prefix already exists, changing prefix to '+new+'\n')
  args.o[0] = new

print("Input files: "+args.f[0]+' '+args.f[1]+' '+args.f[2])
print("Map file: "+args.m[0])
print("Output files: "+args.o[0]+'.ibd12, '+args.o[0]+'.seg')
if args.t[0] == 1:
  print("Printing IBD1 & IBD2 proportions")
else:
  print("Printing IBD1 & IBD2 lengths in cM")
if args.s[0] == 1:
  print("Printing both *.ibd12 and *.seg files\n")
else:
  print("Printing only *.ibd2 file\n")

def mergeIntervals(intervals):
  sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
  merged = []
  
  for higher in sorted_by_lower_bound:
    if not merged:
      merged.append(higher)
    else:
      lower = merged[-1]
      # test for intersection between lower and higher:
      # we know via sorting that lower[0] <= higher[0]
      if higher[0] <= lower[1]:
        upper_bound = max(lower[1], higher[1])
        merged[-1] = (lower[0], upper_bound)  # replace by merged interval
      else:
        merged.append(higher)
  return merged


# Read in map file
loc = {}
for chr in range(1,23):
  loc[chr] = {}

loc_data = open(args.m[0],'r')

for line in loc_data:
  l = str.split(line.rstrip())
  loc[int(l[0])][l[3]] = float(l[2])

loc_data.close()

chr_ends = {}
total_genome = 0
for chr in range(1,23):
  chr_ends[chr] = []
  chr_ends[chr] = [loc[chr][min(loc[chr], key=loc[chr].get)],loc[chr][max(loc[chr], key=loc[chr].get)]]
  total_genome = total_genome + chr_ends[chr][1] - chr_ends[chr][0]


## Collect all segments into single list IBDseg[ind1][ind2][chr]
#run = sys.argv[1]
IBD1_all = {}
IBD2_all = {}
numseg1_all = {}
numseg2_all = {}

for run in range(0,3):
  IBD1 = {}
  IBD2 = {}
  numseg1 = {} #not output
  numseg2 = {} #not output
  IBDseg = {} #pairwise IBD segments: IBDseg[ind1][ind2][chromosome] where value of ind1 < value of ind2
  IBD = open(args.f[run],'r') #open file containing refindIBD output
  for line in IBD:
    l = str.split(line.rstrip())
    #ind1 = str.split(str(l[0]),"_")[0]
    #ind2 = str.split(str(l[2]),"_")[0]
    ind1 = l[0]
    ind2 = l[2]
    h1 = float(l[1])-1 #change homologue values to 0,1 for easier manipulation
    h2 = float(l[3])-1
    if 'chr' in l[4]:
      chr = int(str.split(l[4],'chr')[1])
    else:
      chr = int(l[4])
    if ind1 > ind2: #want ind1 < ind2
      tmp = ind1
      ind1 = ind2
      ind2 = tmp
      tmp = h1
      h1 = h2
      h2 = tmp
    if not ind1 in IBDseg.keys():
      IBDseg[ind1] = {} #create dict for ind1
      IBDseg[ind1][ind2] = {} #create dict for ind2
      IBDseg[ind1][ind2][chr] = [] #create list for chromosome
    elif not ind2 in IBDseg[ind1].keys():
      IBDseg[ind1][ind2] = {} #create dict for ind2
      IBDseg[ind1][ind2][chr] = [] #create list for chromosome
    elif not chr in IBDseg[ind1][ind2].keys():
      IBDseg[ind1][ind2][chr] = [] #create list for chromosome
    IBDseg[ind1][ind2][chr].append([loc[chr][l[5]],loc[chr][l[6]], h1, h2]) #add IBD segment
    #IBDseg[ind1][ind2][chr].append([float(l[5]),float(l[6]), h1, h2]) #add IBD segment
  IBD.close()
  
  for ind1 in IBDseg.keys():
    for ind2 in IBDseg[ind1].keys():
      for chr in IBDseg[ind1][ind2].keys():
        IBDseg[ind1][ind2][chr] = sorted(IBDseg[ind1][ind2][chr]) #sort the segments by start point so we can check for 2 segments overlapping same region (IBD2)
  
  
  IBD1 = {}
  IBD2 = {}
  for ind1 in IBDseg.keys():
    if not ind1 in IBD1.keys():
      IBD1[ind1] = {} #make a list to store each pairwise ind's IBD segs
      IBD2[ind1] = {}
      #numseg1[ind1] = {}
      #numseg2[ind1] = {}
    for ind2 in IBDseg[ind1].keys():
      if not ind2 in IBD1[ind1].keys():
        IBD1[ind1][ind2] = {} #will be the total length of IBD1 segments (in cM) between ind1 and ind2
        IBD2[ind1][ind2] = {} #same as above, but IBD2
        #numseg1[ind1][ind2] = 0 #number of IBD1 segments between ind1 and ind2; not output
        #numseg2[ind1][ind2] = 0 #same as above, but IBD2 segments
      for chr in IBDseg[ind1][ind2].keys():
        if not chr in IBD1[ind1][ind2].keys():
          IBD1[ind1][ind2][chr] = []
          IBD2[ind1][ind2][chr] = []
  
  
  for ind1 in IBDseg.keys():
    for ind2 in IBDseg[ind1].keys():
      for chr in IBDseg[ind1][ind2].keys():
        if len(IBDseg[ind1][ind2][chr]) > 1: #if more than one segment on chromosome 'chr' between ind1 and ind2, check for IBD2
          maxval = len(IBDseg[ind1][ind2][chr])
          k = 0
          while k < maxval - 1:
            kk = k + 1
            while IBDseg[ind1][ind2][chr][k][1] > IBDseg[ind1][ind2][chr][kk][0]: #while first segment ends after a following segment starts
              start = max(IBDseg[ind1][ind2][chr][k][0], IBDseg[ind1][ind2][chr][kk][0])
              end = min(IBDseg[ind1][ind2][chr][k][1], IBDseg[ind1][ind2][chr][kk][1])
              if (abs(IBDseg[ind1][ind2][chr][k][2] - IBDseg[ind1][ind2][chr][kk][2]) == 1 and abs(IBDseg[ind1][ind2][chr][k][3] - IBDseg[ind1][ind2][chr][kk][3]) == 1) or (abs(IBDseg[ind1][ind2][chr][k][2] - IBDseg[ind1][ind2][chr][kk][2]) == 0 and abs(IBDseg[ind1][ind2][chr][k][3] - IBDseg[ind1][ind2][chr][kk][3]) == 0):
                IBD2[ind1][ind2][chr].append([start, end])
              kk = kk + 1
              if kk >= maxval:
                break
            k = k + 1
        IBD2[ind1][ind2][chr] = mergeIntervals(IBD2[ind1][ind2][chr][:])
  
  for ind1 in IBDseg.keys():
    for ind2 in IBDseg[ind1].keys():
      for chr in IBDseg[ind1][ind2].keys():
        IBD1[ind1][ind2][chr] = []
        merged = mergeIntervals(IBDseg[ind1][ind2][chr])
        for i in range(0,len(merged)):
          seg = merged[i]
          k = 0
          if len(IBD2[ind1][ind2][chr]) > 0:
            while  k < len(IBD2[ind1][ind2][chr]):
              while k < len(IBD2[ind1][ind2][chr]) and IBD2[ind1][ind2][chr][k][1] < seg[0]: #while IBD2 is ending before merged segment starts
                k = k + 1
              if k >= len(IBD2[ind1][ind2][chr]):
                IBD1[ind1][ind2][chr].append([seg[0],seg[1]])
                break
              if seg[1] < IBD2[ind1][ind2][chr][k][0]: # if we totally bypassed the IBD1 region
                if k > 0:
                  if IBD2[ind1][ind2][chr][k-1][1] < seg[0]: #if we didn't already add IBD1 from previous IBD2/seg comparison
                    IBD1[ind1][ind2][chr].append([seg[0],seg[1]])
                else:
                  IBD1[ind1][ind2][chr].append([seg[0],seg[1]])
                break
              if k < len(IBD2[ind1][ind2][chr]) - 1:
                if IBD2[ind1][ind2][chr][k][0] >= seg[0] and IBD2[ind1][ind2][chr][k][1] < seg[1] and IBD2[ind1][ind2][chr][k+1][1] < seg[1]: #if this IBD2 and next IBD2 end before merged segment
                  IBD1[ind1][ind2][chr].append([IBD2[ind1][ind2][chr][k][1], IBD2[ind1][ind2][chr][k+1][0]]) #add space between this IBD2 segment and next IBD2 segment
                elif IBD2[ind1][ind2][chr][k][0] >= seg[0] and IBD2[ind1][ind2][chr][k][1] < seg[1] : #if this IBD2 segment ends before merged segment but next IBD2 segment does not
                  #IBD1[ind1][ind2][chr].append([seg[0], IBD2[ind1][ind2][chr][k][0]])
                  IBD1[ind1][ind2][chr].append([IBD2[ind1][ind2][chr][k][1], min(seg[1], IBD2[ind1][ind2][chr][k+1][0])])
              else:
                if IBD2[ind1][ind2][chr][k][0] >= seg[0] and IBD2[ind1][ind2][chr][k][1] < seg[1]: #if this IBD2 segment ends before merged segment but next IBD2 segment does not
                  IBD1[ind1][ind2][chr].append([IBD2[ind1][ind2][chr][k][1], seg[1]])
              if k > 0:
                if IBD2[ind1][ind2][chr][k][0] > seg[0] and IBD2[ind1][ind2][chr][k-1][1] < seg[0]: # if this IBD2 starts after merged segment starts and previous IBD2 segment did not
                  IBD1[ind1][ind2][chr].append([seg[0], IBD2[ind1][ind2][chr][k][0]])
                if seg[1] < IBD2[ind1][ind2][chr][k][0] and IBD2[ind1][ind2][chr][k-1][1] < seg[0]: #if this IBD2 segment starts after merged segment and previous IBD2 segment ends before merged segment
                  IBD1[ind1][ind2][chr].append([seg[0],seg[1]])
              else:
                if IBD2[ind1][ind2][chr][k][0] > seg[0]: # if this IBD2 starts after merged segment starts and there is no previous IBD2 segment
                  IBD1[ind1][ind2][chr].append([seg[0], IBD2[ind1][ind2][chr][k][0]])
              if IBD2[ind1][ind2][chr][k][0] > seg[1]:
                break
              k = k + 1
          else:
            IBD1[ind1][ind2][chr].append([seg[0],seg[1]]) 
  
  for ind1 in IBD1.keys():
    for ind2 in IBD1[ind1].keys():
      for chr in IBD1[ind1][ind2].keys():
        if not ind1 in IBD1_all.keys():
          IBD1_all[ind1] = {}
        if not ind2 in IBD1_all[ind1].keys():
          IBD1_all[ind1][ind2] = {}
        if not chr in IBD1_all[ind1][ind2].keys():
          IBD1_all[ind1][ind2][chr] = []
        if len(IBD1[ind1][ind2][chr]) > 0:
          for k in range(0, len(IBD1[ind1][ind2][chr])):
            IBD1_all[ind1][ind2][chr].append(IBD1[ind1][ind2][chr][k])
  
  for ind1 in IBD2.keys():
    for ind2 in IBD2[ind1].keys():
      for chr in IBD2[ind1][ind2].keys():
        if not ind1 in IBD2_all.keys():
          IBD2_all[ind1] = {}
        if not ind2 in IBD2_all[ind1].keys():
          IBD2_all[ind1][ind2] = {}
        if not chr in IBD2_all[ind1][ind2].keys():
          IBD2_all[ind1][ind2][chr] = []
        if len(IBD2[ind1][ind2][chr]) > 0:
          for k in range(0, len(IBD2[ind1][ind2][chr])):
            IBD2_all[ind1][ind2][chr].append(IBD2[ind1][ind2][chr][k])


for ind1 in IBD1_all.keys():
  for ind2 in IBD1_all[ind1].keys():
    for chr in IBD1_all[ind1][ind2].keys():
      IBD1_all[ind1][ind2][chr] = mergeIntervals(IBD1_all[ind1][ind2][chr][:])

for ind1 in IBD2_all.keys():
  for ind2 in IBD2_all[ind1].keys():
    for chr in IBD2_all[ind1][ind2].keys():
      IBD2_all[ind1][ind2][chr] = mergeIntervals(IBD2_all[ind1][ind2][chr][:])


#IBD2_new = IBD2_all.copy()
IBD1_new = {}
for ind1 in IBD1_all.keys():
  if not ind1 in IBD1_new.keys():
    IBD1_new[ind1] = {}
  for ind2 in IBD1_all[ind1].keys():
    if not ind2 in IBD1_new[ind1].keys():
      IBD1_new[ind1][ind2] = {}
    for chr in IBD1_all[ind1][ind2].keys():
      if not chr in IBD1_new[ind1][ind2].keys():
        IBD1_new[ind1][ind2][chr] = []
      k = 0 #IBD1
      kk = 0  #IBD2
      merged = mergeIntervals(IBD1_all[ind1][ind2][chr])
      if len(IBD2_all[ind1][ind2][chr]):
        while kk < len(IBD2_all[ind1][ind2][chr]) and k < len(merged):
          while k < len(merged) and merged[k][1] < IBD2_all[ind1][ind2][chr][kk][0]: #while our IBD1 segment ends before IBD2 segment starts
            if not [merged[k][0],merged[k][1]] in IBD1_new[ind1][ind2][chr]:
              if len(IBD1_new[ind1][ind2][chr]):
                if merged[k][0] > IBD1_new[ind1][ind2][chr][-1][1]:
                  IBD1_new[ind1][ind2][chr].append(merged[k])
              else:
                IBD1_new[ind1][ind2][chr].append(merged[k])
            k = k + 1
          if k < len(merged):
            while kk < len(IBD2_all[ind1][ind2][chr]) and IBD2_all[ind1][ind2][chr][kk][1] <= merged[k][0]: #while our IBD2 segments ends before IBD1 segment starts
              kk = kk + 1
            if kk < len(IBD2_all[ind1][ind2][chr]):
              if merged[k][0] == IBD2_all[ind1][ind2][chr][kk][0] and merged[k][1] == IBD2_all[ind1][ind2][chr][kk][1]: #if the segments are identical
                kk = kk + 1 #move to next IBD2 segment
              elif merged[k][1] <= IBD2_all[ind1][ind2][chr][kk][0]: #end of IBD1 segment occurs before IBD2 segment
                if len(IBD1_new[ind1][ind2][chr]):
                  if merged[k][0] > IBD1_new[ind1][ind2][chr][-1][1]:
                    IBD1_new[ind1][ind2][chr].append(merged[k])
                else:
                  IBD1_new[ind1][ind2][chr].append(merged[k])
                k = k + 1
              elif merged[k][0] <= IBD2_all[ind1][ind2][chr][kk][0] and merged[k][1] >= IBD2_all[ind1][ind2][chr][kk][1]: #if the IBD2 segment is contained within IBD1 segment
                if len(IBD1_new[ind1][ind2][chr]):
                  start = max(merged[k][0],IBD1_new[ind1][ind2][chr][-1][1])
                  if start != IBD2_all[ind1][ind2][chr][kk][0]:
                    IBD1_new[ind1][ind2][chr].append([start, IBD2_all[ind1][ind2][chr][kk][0]])
                else:
                  IBD1_new[ind1][ind2][chr].append([merged[k][0],IBD2_all[ind1][ind2][chr][kk][0]])
                if kk < len(IBD2_all[ind1][ind2][chr]) - 1:
                  end = min(merged[k][1],IBD2_all[ind1][ind2][chr][kk+1][0])
                  IBD1_new[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][kk][1],end])
                  if end == merged[k][1]:
                    k = k + 1
                else:
                  IBD1_new[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][kk][1],merged[k][1]])
                kk = kk + 1
              elif merged[k][0] <= IBD2_all[ind1][ind2][chr][kk][0] and merged[k][1] <= IBD2_all[ind1][ind2][chr][kk][1]: #if the IBD2 segment starts after IBD1 segment and IBD2 segment ends after IBD1 segment
                #IBD1_new[ind1][ind2][chr].append([merged[k][0],IBD2_all[ind1][ind2][chr][kk][0]]) #old
                if len(IBD1_new[ind1][ind2][chr]):
                  IBD1_new[ind1][ind2][chr].append([max(merged[k][0], IBD1_new[ind1][ind2][chr][-1][1]), min(IBD2_all[ind1][ind2][chr][kk][0],merged[k][1])])
                else:
                  IBD1_new[ind1][ind2][chr].append([merged[k][0], min(IBD2_all[ind1][ind2][chr][kk][0], merged[k][1])])
                k = k + 1
              elif merged[k][0] >= IBD2_all[ind1][ind2][chr][kk][0] and merged[k][1] >= IBD2_all[ind1][ind2][chr][kk][1]: #if the IBD2 segment starts before IBD1 segment and IBD2 segment ends before IBD1 segment
                if kk < len(IBD2_all[ind1][ind2][chr]) - 1:
                  IBD1_new[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][kk][1],min(merged[k][1],IBD2_all[ind1][ind2][chr][kk+1][0])])
                else:
                  IBD1_new[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][kk][1],merged[k][1]])
                kk = kk + 1
              else: #if the IBD2 segment contains IBD1 segment, we don't care
                k = k + 1
        
        while k < len(merged):
          if len(IBD1_new[ind1][ind2][chr]):
            if IBD1_new[ind1][ind2][chr][-1][1] <= merged[k][0]:
              IBD1_new[ind1][ind2][chr].append(merged[k])
            k = k + 1
          else:
            IBD1_new[ind1][ind2][chr].append(merged[k])
            k = k + 1
      else:
        for k in range(0,len(IBD1_all[ind1][ind2][chr])):
          IBD1_new[ind1][ind2][chr].append(IBD1_all[ind1][ind2][chr][k])
        
        


IBD1 = {}
IBD2 = {}
for ind1 in IBD1_new.keys():
  for ind2 in IBD1_new[ind1].keys():
    if not ind1 in IBD1.keys():
      IBD1[ind1] = {}
      IBD2[ind1] = {}
    if not ind2 in IBD1[ind1].keys():
      IBD1[ind1][ind2] = 0
      IBD2[ind1][ind2] = 0
    for chr in IBD1_new[ind1][ind2].keys():
      if len(IBD1_new[ind1][ind2][chr]) > 0:
        for k in range(0, len(IBD1_new[ind1][ind2][chr])):
          IBD1[ind1][ind2] = IBD1[ind1][ind2] + (IBD1_new[ind1][ind2][chr][k][1] - IBD1_new[ind1][ind2][chr][k][0])

for ind1 in IBD2_all.keys():
  for ind2 in IBD2_all[ind1].keys():
    if not ind1 in IBD2.keys():
      IBD1[ind1] = {}
      IBD2[ind1] = {}
    if not ind2 in IBD2[ind1].keys():
      IBD1[ind1][ind2] = 0
      IBD2[ind1][ind2] = 0
    for chr in IBD2_all[ind1][ind2].keys():
      if len(IBD2_all[ind1][ind2][chr]) > 0:
        for k in range(0, len(IBD2_all[ind1][ind2][chr])):
          IBD2[ind1][ind2] = IBD2[ind1][ind2] + (IBD2_all[ind1][ind2][chr][k][1] - IBD2_all[ind1][ind2][chr][k][0])


if args.t[0] == 1:
  for ind1 in IBD1.keys():
    for ind2 in IBD1[ind1].keys():
      IBD1[ind1][ind2] = IBD1[ind1][ind2]/total_genome
      IBD2[ind1][ind2] = IBD2[ind1][ind2]/total_genome

  
outfile_IBD = open(args.o[0]+'.ibd12','w')
for ind1 in IBD1.keys():
  for ind2 in IBD1[ind1].keys():
    outfile_IBD.write(ind1+'\t'+ind2+'\t'+str(IBD1[ind1][ind2])+'\t'+str(IBD2[ind1][ind2])+'\n') #outputs ind1 ind2 IBD1_between_ind1_ind2 IBD2_between_ind1_ind2

outfile_IBD.close()

if args.s[0] == 1:
  outfile_IBD = open(args.o[0]+'.seg','w')
  for ind1 in IBD1_new.keys():
    for ind2 in IBD1_new[ind1].keys():
      for chr in IBD1_new[ind1][ind2].keys():
        for i in range(0,len(IBD1_new[ind1][ind2][chr])):
          outfile_IBD.write(ind1+'\t'+ind2+'\t'+str(chr)+'\tIBD1\t'+str(IBD1_new[ind1][ind2][chr][i][0])+'\t'+str(IBD1_new[ind1][ind2][chr][i][1])+'\n')

  for ind1 in IBD2_all.keys():
    for ind2 in IBD2_all[ind1].keys():
      for chr in IBD2_all[ind1][ind2].keys():
        for i in range(0, len(IBD2_all[ind1][ind2][chr])):
          outfile_IBD.write(ind1+'\t'+ind2+'\t'+str(chr)+'\tIBD2\t'+str(IBD2_all[ind1][ind2][chr][i][0])+'\t'+str(IBD2_all[ind1][ind2][chr][i][1])+'\n')

  outfile_IBD.close()



