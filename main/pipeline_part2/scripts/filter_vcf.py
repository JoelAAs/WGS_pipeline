import pandas as pd
import argparse
import  sys

class FilterVCF:
    def __init__(self, vcf_file, snp_list, sample_file, output, chunksize=100000, exclude_snp=True,
                 exclude_samples=True):
        self.vcffile = open(vcf_file, "r")
        self.snps = self.read_snp_list(snp_list)
        self.samples = (self.read_sample_list(sample_file) if exclude_samples else [])
        self.chunksize = chunksize
        self.exclude_snp = exclude_snp
        # self.exclude_samples = exclude_samples
        self.variant_header = None
        self.output = open(output, "w")
        self.id_col = ["CHROM", "POS", "REF", "ALT"]
        self.get_chrom_order = {
            "1": 1,  "2": 2,  "3": 3,  "4": 4,  "5": 5,  "6": 6,  "7": 7,  "8": 8,  "9": 9,  "10": 10,
            "11": 11, "12": 12, "13": 13, "14": 14, "15": 15, "16": 16, "17": 17, "18": 18, "19": 19, "20": 20,
            "21": 21, "22": 22, "X": 23, "Y": 24, "MT": 25, "GL000207.1": 26, "GL000226.1": 27, "GL000229.1": 28,
            "GL000231.1": 29, "GL000210.1": 30, "GL000239.1": 31, "GL000235.1": 32, "GL000201.1": 33, "GL000247.1": 34,
            "GL000245.1": 35, "GL000197.1": 36, "GL000203.1": 37, "GL000246.1": 38, "GL000249.1": 39, "GL000196.1": 40,
            "GL000248.1": 41, "GL000244.1": 42, "GL000238.1": 43, "GL000202.1": 44, "GL000234.1": 45, "GL000232.1": 46,
            "GL000206.1": 47, "GL000240.1": 48, "GL000236.1": 49, "GL000241.1": 50, "GL000243.1": 51, "GL000242.1": 52,
            "GL000230.1": 53, "GL000237.1": 54, "GL000233.1": 55, "GL000204.1": 56, "GL000198.1": 57, "GL000208.1": 58,
            "GL000191.1": 59, "GL000227.1": 60, "GL000228.1": 61, "GL000214.1": 62, "GL000221.1": 63, "GL000209.1": 64,
            "GL000218.1": 65, "GL000220.1": 66, "GL000213.1": 67, "GL000211.1": 68, "GL000199.1": 69, "GL000217.1": 70,
            "GL000216.1": 71, "GL000215.1": 72, "GL000205.1": 73, "GL000219.1": 74, "GL000224.1": 75, "GL000223.1": 76,
            "GL000195.1": 77, "GL000212.1": 78, "GL000222.1": 79, "GL000200.1": 80, "GL000193.1": 81, "GL000194.1": 82,
            "GL000225.1": 83, "GL000192.1": 84
        }
        self.eof = False

    def read_snp_list(self, snp_list):
        with open(snp_list, "r") as f:
            snps = [l.strip().split(":") for l in f.readlines()]

        snp_df = pd.DataFrame(snps, columns=self.id_col)
        snp_df["chrom_order"] = snp_df["CHROM"].apply(lambda x: self.get_chrom_order[x], dtype={"POS":int})
        snp_df = snp_df.sort_values(["POS", "chrom_order"])
        snp_df["snps_id"] = snp_df[self.id_col].apply(lambda x: ":".join(x.values.astype(str)))

        return list(snp_df["snps_id"])

    def read_sample_list(self, sample_file):
        with open(sample_file, "r") as f:
            samples = [line.strip() for line in f.readlines()]

        return samples

    def set_meta_and_get_header(self):
        for line in self.vcffile:
            if "#CHROM" in line:
                self.variant_header = line.strip()[1:].split("\t")
                print(self.variant_header)
                vcf_header = [head for head in self.variant_header if head not in self.samples]
                self.output.write("#" + "\t".join(vcf_header) + "\n")
                break
            else:
                self.output.write(line)

    def _read_chunk(self, chunksize):
        lines = [""] * chunksize
        line_count = 0
        while line_count < chunksize:
            try:
                lines[line_count] = next(self.vcffile)
                line_count += 1
            except StopIteration:
                self.eof = True
                lines = filter(lambda x: x != "", lines)
                break

        lines = list(map(lambda x: x.split("\t"), lines))

        return pd.DataFrame(lines, columns=self.variant_header)

    def filter(self):
        self.set_meta_and_get_header()
        while not self.eof:
            current_chunk_df = self._read_chunk(self.chunksize)
            if self.samples:
                for sample in self.samples:
                    try:
                        del current_chunk_df[sample]
                    except KeyError as e:
                        raise e

            if self.snps:
                del_snp = False
                chunk_iter = current_chunk_df.iterrows()
                filter_chrom, filter_pos, filter_ref, filter_alt = self.snps[0].split(":")
                filter_chrom_order = self.get_chrom_order[filter_chrom]

                while True:
                    try:
                        _, row = next(chunk_iter)
                    except StopIteration:
                        break
                    if self.snps:
                        if del_snp:
                            print(self.snps[0].split(":"))
                            filter_chrom, filter_pos, filter_ref, filter_alt = self.snps[0].split(":")
                            filter_chrom_order = self.get_chrom_order[filter_chrom]
                            del_snp = False

                        print(row[self.id_col])
                        c_chrom, c_pos, c_ref, c_alt = row[self.id_col]
                        c_chrom_order = self.get_chrom_order[c_chrom]

                        if filter_chrom_order > c_chrom_order or \
                                filter_pos > c_pos or \
                                c_ref != filter_ref or \
                                c_alt != filter_alt:  # We assume that all variants in list are present in file
                            self.output.write("\t".join(row.values.astype(str)) + "\n")
                        else:
                            del self.snps[0]
                            del_snp = True
                    else:
                        self.output.write("\t".join(row.values.astype(str)) + "\n")
            else:
                chunk_iter = current_chunk_df.iterrows()
                for _, row in chunk_iter:
                    self.output.write("\t".join(row.values.astype(str)) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Remove SNPs and samples from vcf"
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF to filter out variants from"
    )
    parser.add_argument(
        "--snpslist",
        type=str,
        required=True,
        help="List without header containing variants on form 'CHROM:POS:REF:ALT"
    )
    parser.add_argument(
        "--samplelist",
        type=str,
        required=True,
        help="List without header containing samples to be excluded"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to new filtered VCF"
    )

    args = parser.parse_args(sys.argv[1:])

    filt = FilterVCF(
        args.vcf,
        args.snplist,
        args.samplelist,
        args.output,
        chunksize=120000  # should be ~58.6 GB of variants
    )

    filt.filter()


if __name__ == '__main__':
    main()
