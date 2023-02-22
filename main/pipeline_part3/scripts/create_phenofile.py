import sys
import re
import argparse
import collections


def mkParser():
  parser = argparse.ArgumentParser(description = "Create phenofile for plink")
  parser.add_argument("--samplefile",   type = str,    required = True,   help="the main samplefile")
  parser.add_argument("--out",          type = str,    required = True,   help="the name of the output file")

  return parser.parse_args()
			

# create sets with which adrs and drug types present in file
def make_dictionary(infile):
  fh         = open(infile, mode = 'r')
  line       = fh.readline()
  sampledict = {}
  adrset     = set()
  drugset    = set()

  for line in fh:
    cols = line.strip().split()
    sampledict[cols[0]] = {'adrtype': [], 'drugtype': []}

    if len(cols) > 6: # watch this lenght might change in file
      adr  = cols[5].strip().split(';')  
      drug = cols[6].strip().split(';')  # this needs to be changed later

      [sampledict[cols[0]]['adrtype'].append(i) for i in adr]
      [sampledict[cols[0]]['drugtype'].append(j) for j in drug]
     
      if adr != "":
        sampledict[cols[0]]['drugtype'].append("Case") # Sorry for adding ugly code /Joel

      [adrset.add(i) for i in adr]
      [drugset.add(j) for j in drug]
  
  drugset.add("Case")
  return sampledict, adrset, drugset


# create all possible combinations for assoc testing
def make_comparisons(adrset, drugset):
  assoc_list = []
  for adr in adrset:
    assoc_list.append(adr+'~rest')
    assoc_list.append(adr+'~control')
  for drug in drugset:
    assoc_list.append(drug+'~control')
    assoc_list.append(drug+'~rest')

  for drug in drugset:
    for adr in adrset:
      assoc_list.append(drug+'_'+adr+'~'+drug+'_rest')

  return assoc_list


# write case control status for each test
def identify(sample, assoc_list):
  status_list = []
  for assoc in assoc_list:
    cols = re.split('_|~',assoc)

    #check for adr group and drug comparisons only
    if len(cols) == 2:
      if cols[0] in sample['adrtype'] or cols[0] in sample['drugtype']:
        status_list.append(2)
      elif cols[1] == 'control' and len(sample['adrtype']) > 0:
        status_list.append(-9)
      else:
        status_list.append(1)

    #check for combinations of drug and adr
    else:
      if cols[1] in sample['adrtype'] and cols[0] in sample['drugtype']:
        status_list.append(2)
      elif cols[0] in sample['drugtype']:
        status_list.append(1)
      else:
        status_list.append(-9)

  return status_list


# create final dict for writing to file
def make_assocdict(sdict, assoc_list):
  count_dict = {}
  for sample in sdict:
    count_dict[sample] = identify(sdict[sample], assoc_list)

  return count_dict


# summarize counts for each assoc analysis and write to file 
def write_counttable(count_dict, assoc_list, out):
  summary_dict = {}
  i = 0
  while i < len(assoc_list):
    summary_dict[assoc_list[i]] = {'case': 0, 'control':0, 'missing': 0}
    for sample in count_dict:
      if count_dict[sample][i] == 2:
        summary_dict[assoc_list[i]]['case'] += 1
      elif count_dict[sample][i] == 1:
        summary_dict[assoc_list[i]]['control'] += 1
      elif count_dict[sample][i] == -9:
        summary_dict[assoc_list[i]]['missing'] += 1
    i += 1
  sort = sorted(summary_dict.keys(), key=lambda x: (summary_dict[x]['case'], summary_dict[x]['control']), reverse=True)
  index_map = {v: i for i, v in enumerate(sort)}
  
  out = open(out+'_counts.table', 'w')
  out.write('test\tcases\tcontrols\tmissing\n')
  for total in sorted(summary_dict.items(), key=lambda pair: index_map[pair[0]]):
    out.write("{}\t{}\t{}\t{}\n".format(total[0],total[1]['case'],total[1]['control'],total[1]['missing']))
  out.close()


def write_outfile(assoc_list, count_dict, out):
  out = open(out, 'w')
  out.write('IID')
  for s in assoc_list:
    out.write('\t'+s)
  out.write('\n')

  for sample in sorted(count_dict):
    out.write(sample)
    for ass in count_dict[sample]:
      out.write('\t'+str(ass))
    out.write('\n')
  out.close()


def main():
  print("### INFO ###   Running...")
  args                   = mkParser()  
  sdict, adrset, drugset = make_dictionary(args.samplefile)
  assoc_list             = make_comparisons(adrset, drugset)
  count_dict             = make_assocdict(sdict, assoc_list)
  write_counttable(count_dict, assoc_list, args.out)
  write_outfile(assoc_list, count_dict, args.out)
  print("### INFO ###   Done!")

main()
