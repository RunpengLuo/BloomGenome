import sys
import os

from utils.vcf_reader import VCF_Reader
from utils.kinship_inference import kinship_inference_terry, gtruth
from utils.bloom_filter import *

if __name__ == "__main__":
    # print("{0} <vcf_file> <outdir> <samples> <overwrite_mode>")
    # [_, vcf, samplef, outdir] = sys.argv
    # # overwrite = sys.argv[4] == '1' if len(sys.argv) == 5 else False

    # if outdir[-1] != "/":
    #     outdir += "/"
    
    # if not os.path.exists(outdir):
    #     os.system("mkdir " + outdir)
    #     print("output directory does not exist, creating..")
    # else:
    #     print("output directory exists, using intermediate computed results")

    # vcf_reader = VCF_Reader(vcf, samplef, outdir)
    # vcf_reader.split_samples()
    # print("Time for step 1: ", vcf_reader.step1time)
    # print("Time for step 2: ", vcf_reader.step2time)



    kinship_inference_terry(
        "chr22test_2011/NA19660_0.txt", "chr22test_2011/NA19660_1.txt", "chr22test_2011/NA19660_pfilter.txt",
        "chr22test_2011/NA19685_0.txt", "chr22test_2011/NA19685_1.txt", "chr22test_2011/NA19685_pfilter.txt",
        1000000
        )

    gtruth(
        "chr22test_2011/NA19660_0.txt", "chr22test_2011/NA19660_1.txt", "chr22test_2011/NA19660_pfilter.txt",
        "chr22test_2011/NA19685_0.txt", "chr22test_2011/NA19685_1.txt", "chr22test_2011/NA19685_pfilter.txt")

    




