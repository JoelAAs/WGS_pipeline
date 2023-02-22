from scipy.stats   import uniform, randint, chi2
from statistics    import median
from os            import path

import sys
sys.path.insert(0, "/proj/sens2018106/softwares/python_packages")

import argparse
import numpy             as np
import pandas            as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def mkParser():
	parser     = argparse.ArgumentParser(description = "Makes Manhattan and QQ plots from plink assoc files")
	parser.add_argument("--assoc",        type = str,    required = True,    help="output from plink assoc or logistic")
	parser.add_argument("--out",          type = str,    required = True,    help="the name of the output file")

	return parser.parse_args()
    
    
    
def readfile(assoc):
    fh      = open(assoc, "r")
    line    = fh.readline()
    header  = line.strip().split()
    chr_ix  = header.index('#CHROM')
    bp_ix   = header.index('POS')
    p_ix    = header.index('P')
    test_ix = header.index('TEST')

    s_dict = {}
    ids    = 1
    allowed_chroms = list(map(str, range(1,23)))
    allowed_chroms.append("X")
    for line in fh:
        cols = line.strip().split()
        if cols[test_ix] == 'ADD' and cols[p_ix] != 'NA' and cols[chr_ix] in allowed_chroms:     # only select test ADD and rows with non-empty pvalues
            s_dict[ids] = {'#CHROM': (int(cols[chr_ix]) if cols[chr_ix] != "X" else 23), 
                           'POS': int(cols[bp_ix]),
                           'P': float(cols[p_ix]), 
                           '-log(P)': -np.log10(float(cols[12])),  # idx 11 for not firth
                           'chisq': chi2.ppf(1-float(cols[12]), df=1)}  # this one takes long
            ids += 1
    df = pd.DataFrame.from_dict(s_dict, orient = 'index')
    print(df.columns.values) 
    df["#CHROM"] = df["#CHROM"].astype('category')
    df = df.sort_values('#CHROM')
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('#CHROM'))
    
    return df, df_grouped


def plot(df, df_grouped, filename, out):
    test = ''.join(path.basename(filename).split('.')[0])
    grid = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3)

    plt.figure(figsize=(20,20))
    ax = plt.subplot(grid[0,:2])
    colors = ['blue','grey']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='-log(P)',color=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, size = 15)
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, df['-log(P)'].max()+1])
    ax.set_xlabel('Chromosome', size = 20)
    ax.set_ylabel('-log10(P)', size = 20)
    plt.axhline(y=7.3, color='r', linestyle='-')
    plt.axhline(y=5, color='b', linestyle='-')
    plt.title(test, size = 30)

    lambdaGC = median(df.chisq)/chi2.ppf(0.5, df=1)

    plt.subplot(grid[1,0])
    observed = sorted(df['P'])
    lobs     = [-np.log10(i) for i in observed]

    expected = list(range(1,len(observed)+1))
    lexp = [-np.log10((i) / (len(expected)+1)) for i in expected]

    plt.plot(lexp, lobs, marker='o', markersize = 2, color = 'b', linestyle=' ')
    plt.plot(list(range(0,int(round(max(lexp)))+1)), list(range(0,int(round(max(lexp)))+1)), color = 'darkred')
    plt.text(max(lexp)-1.5, 1, 'LambdaGC: \n'+str(round(lambdaGC,2)),
            verticalalignment='top', size = 20)
    plt.title('QQ-plot', size = 20)
    plt.xlabel('Expected -log10(P)')
    plt.ylabel('Observed -log10(P)')

    plt.savefig(out)
    plt.close()



def main():
    print('### INFO ###  Plotting...')
    args = mkParser()
    df, df_grouped = readfile(args.assoc)
    plot(df, df_grouped, args.assoc, args.out)
    print('### INFO ###  Done!')


if __name__ == "__main__":
	main()
