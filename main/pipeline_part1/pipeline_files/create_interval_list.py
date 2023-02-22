import sys, os, glob
import argparse

""" 
Utility script to create the list of section that the genome will be divided in for parallel processing
"""

def main(args):
    if args.chr and args.block is not None:
        sys.exit("ERROR: only one between block and chr can be specified")
    #memorise the dict, assume it ordered
    chr_list = parse_dict(args.dict, args.outfolder)
    if args.chr:
        chr_intervals(chr_list, args.outfolder)
    else:
        intervals = block_intervals(chr_list, args.block, args.outfolder)



def block_intervals(chr_list, block_size, outfolder):
    blocks = []
    block  = []
    current_block_left     = block_size # initialise this
    for chr in chr_list:
        current_chr        = chr[0]
        current_chr_length = chr[1]
        current_interval_start = 1 # beginning of a chr, need to reset the interval start
        while current_interval_start + current_block_left <= current_chr_length:
            block.append([current_chr,current_interval_start, current_interval_start + current_block_left]) # this block if full
            blocks.append(block)
            current_interval_start = current_interval_start + current_block_left # next block starts right after the last
            block = [] # reset block
            current_block_left = block_size # all the block is left
        # need to handle the  last segment of the chr or maybe the all chr if it is particularly small
        if current_interval_start < current_chr_length: #check special case when the last inserted block overalps the exact last base of a chr
            block.append([current_chr,current_interval_start, current_chr_length]) #fill in the block with this interval
            blocks.append(block)
            block = []
            current_block_left = block_size
    #now I need to add the last block that might have been not added
    if len(block) > 0:
        blocks.append(block)
    #now proint the intervals
    cwd      = os.getcwd()
    current_interval = 1

    with open(outfolder + "00_intervals.txt", "w") as intervals_file:
        for block in blocks:
            interval_file_name = os.path.join(cwd, outfolder + 'interval{0:04d}.intervals'.format(current_interval))
            current_interval += 1
            intervals_file.write("{}\n".format(interval_file_name))
            with open(interval_file_name, "w") as interval_file:
                for section in block:
                    interval_file.write("{}:{}-{}\n".format(section[0],section[1],section[2]))




def chr_intervals(chr_list, outfolder):
    MT_GL = []
    intervals = []
    for chr in chr_list:
        if "GL" not in chr[0] and "MT" not in chr[0]:
            intervals.append([chr[0]])
        else:
            MT_GL.append(chr[0])
    if len(MT_GL) > 0:
        intervals.append(MT_GL)
    cwd      = os.getcwd()
    current_interval = 1
    with open(outfolder + "00_intervals.txt", "w") as intervals_file:
        for interval in intervals:
            interval_file_name = os.path.join(cwd, outfolder + 'interval{0:04d}.intervals'.format(current_interval))
            current_interval += 1
            intervals_file.write("{}\n".format(interval_file_name))
            with open(interval_file_name, "w") as interval_file:
                for section in interval:
                    interval_file.write("{}\n".format(section))


        



def parse_dict(chr_sam_file, outfolder):
    if not os.path.exists(chr_sam_file):
        sys.exit("ERROR: dict file provided does not exists.")
    chr_list = []
    with open(chr_sam_file, "r") as chr_sam:
        for sam_line in chr_sam:
            sam_list = sam_line.split("\t")
            if "@SQ" not in sam_list[0]:
                continue
            chr = sam_list[1].split(":")[1]
            chr_length = int(sam_list[2].split(":")[1])
            chr_list.append([chr, chr_length])
    return chr_list



if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Utility script to create the list of section that the genome will be divided in for parallel processing. The dict is assumed ordered, i.e., chr1m chr2, etc.""")
    parser.add_argument('--dict',      type = str,   help = "dictionary of the genome (in general a .dict file) in SAM fomat (i.e. @SQ SN:NAME LN:LENGTH UR:OTHER_STUFF", required=True)
    parser.add_argument('--block',     type = int,   help = "Size of the blocks the genome will be divided up")
    parser.add_argument('--chr',       help = "each chr is a block (assume human genome, so threat as single block unk and alternate", action='store_true', default=False)
    parser.add_argument('--outfolder', type = str,   help = "The name of the output folder in which to write the results")
    args = parser.parse_args()
    main(args)
