import sys
import re
import argparse


def mkParser():
  parser = argparse.ArgumentParser(description = "Translates the sequencing output sample name to real sample names and gender using a sample file")
  parser.add_argument("--bfile",        type = str,    required = True,   help="the plink file to convert")
  parser.add_argument("--sample",       type = str,    required = True,   help="a samplefile with gender information for each sample. Sample name in column 2 and gender in column 3")
  parser.add_argument("--out",          type = str,    required = True,   help="the name of the output file")


  return parser.parse_args()
			
def make_dictionary(infile):
  sampleDict = {}
  for line in infile:
    cols = line.strip().split('\t')
    sampleDict[cols[1]] = cols[2]
  print(sampleDict)
  return sampleDict

def translate_samples(fam, dic, out):
  out = open(out+".fam", "w")
  sep = ' '
  for line in fam:
    cols = line.strip().split(' ')
    if not cols[1] in dic:
      print(str(cols[1])+'  not in sampledict. Exiting!')
      exit()
    out.write(str(cols[0])+sep+str(cols[1])+sep+str(cols[2])+sep+str(cols[3])+sep+str(dic[cols[1]])+sep+str(cols[5])+'\n')
  out.close()      

def translate_variants(bim, out):
  out = open(out+".bim", "w")
  sep = ' '
  for line in bim:
    cols = line.strip().split()
    var  = str(cols[0])+'_'+str(cols[3])+'_'+str(cols[4])+':'+str(cols[5])
    out.write(str(cols[0])+sep+var+sep+sep.join(cols[2:])+'\n')
  out.close()

def main():
  args = mkParser()  
  print("##  INFO  ###   Writing new fam file")
  infile = open(args.sample, "r")
  dic    = make_dictionary(infile)
  fam    = open(args.bfile+".fam", "r")
  bim    = open(args.bfile+".bim","r")
  translate_samples(fam, dic, args.out)
  translate_variants(bim, args.out)
  print("##  info  ###   Done!")

main()
