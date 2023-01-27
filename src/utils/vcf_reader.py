import os
import re
import sys
from bitarray import bitarray
import time

CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = range(9)
R, A1, A2 = range(3)
UNPHASED, PHASED = '/', '|' 

# Original Code: Terry Yu
# Reformat and Optimize using bcftools: John Luo

class VCF_Reader:
    def __init__(self, vcf_file: str, samples_file: str, outdir: str) -> None:
        print("Extracting sample ids from " + samples_file)
        t1s = time.time()

        sfd = open(samples_file, "r")
        samples = list(set(sid.strip() for sid in sfd.readlines()))
        sfd.close()

        succ_samples = []

        print("Extract samples: {0} from vcf file".format(samples))
        
        bcf_file = outdir + "data.bt.bcf"
        if os.path.exists(bcf_file):
            print("BCF file exists")
        else:
            print("Convert VCF file to filtered compressed BCF file..")
            os.system("bcftools view --force-samples -s {0} -O b {1} > {2}".format(','.join(samples), vcf_file, bcf_file))
        
        pra_file = outdir + "pra.txt"
        if os.path.exists(pra_file):
            print("POS+REF+ALT file exists")
        else:
            print("Obtain POS+REF+ALT information..")
            os.system("bcftools query -f '%POS\t%REF\t%ALT\n' {0} -o {1}".format(bcf_file, pra_file))

        rcode = 0
        for sample in samples:
            sample_file = "{0}{1}.txt".format(outdir, sample)
            print("Extracting information for sample: " + sample)
            rcode = os.system("bcftools view -O b -s {0} {1} | bcftools query -f '[%GT]\n' -o {2}".format(
                sample, bcf_file, sample_file
            ))
            if rcode != 0:
                print("Fail to extract sample: " + sample)
            else:
                print("Success to extract sample: " + sample)
                succ_samples.append(sample)
                
                rfd = open(sample_file, "r")
                # assume alt index < 10
                pfilter = bitarray(map(lambda x: x[0] != '.' and x[2] != '.', rfd.readlines()))
                rfd.close()

                pfilter_file = "{0}{1}_pfilter.txt".format(outdir, sample)
                os.system("echo "" > " + pfilter_file)                
                
                pfd = open(pfilter_file, "w")
                pfd.write("{0}\n".format(pfilter.to01()))
                pfd.close()

        self.succ_samples = succ_samples
        self.outdir = outdir
        self.step1time = round(time.time() - t1s, 2)
        self.step2time = 0

        print("VCF reader initialized")
        print("All the successful extracted samples: ", succ_samples)
        return
    
    def split_samples(self):
        tss = time.time()
        pra_list = []
        pra_file = self.outdir + "pra.txt"
        max_variation = 0
        with open(pra_file, "r") as pfd:
            for line in pfd:
                pos, ref, alts = line.strip().split("\t")
                alt_list = alts.split(",")
                
                max_variation = max(max_variation, len(alt_list))
                
                alt_list.extend([ref, '.'])
                pra_list.append((pos, ref, alt_list))
            pfd.close()

        # for fast lookup instead of type conversion
        roundup_idx = {}
        for i in range(1, max(2, max_variation)):
            roundup_idx[str(i)] = i - 1
        # special case for ref
        roundup_idx['0'] = -2
        # special case for missing
        roundup_idx['.'] = -1

        for sample in self.succ_samples:
            ts = time.time()
            sample_file = "{0}{1}.txt".format(self.outdir, sample)
            afile0 = "{0}{1}_0.txt".format(self.outdir, sample)
            afile1 = "{0}{1}_1.txt".format(self.outdir, sample)
            os.system("echo "" > {0}; echo "" > {1}".format(afile0, afile1))
            afd1, afd2 = open(afile0, "w"), open(afile1, "w")
            sfd = open(sample_file, "r")
            for (pos, ref, alts), line in zip(pra_list, sfd.readlines()):
                try:
                    a1, a2 = re.split("[\|]", line.strip())
                    afd1.write("{0}:{1}:{2}\n".format(pos, ref, alts[roundup_idx[a1]]))
                    afd2.write("{0}:{1}:{2}\n".format(pos, ref, alts[roundup_idx[a2]]))
                except:
                    print(roundup_idx)
                    sys.exit(1)

            afd1.close()
            afd2.close()
            sfd.close()
            te = time.time()
            print("Sample: {0} finished, time elapsed: {1}".format(sample, round(te-ts, 2)))
        
        self.step2time = round(time.time() - tss, 2)
        return 0
            