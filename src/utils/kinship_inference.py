from utils.bloom_filter import *
from utils.utils import *
from bitarray import bitarray

import sys

relationship_dict = {
     'monozygotic_twin': (sys.maxsize, 0.35355339059327373) # > 1/(2**1.5)
    ,'parent_offspring_or_full_sib': (0.35355339059327373, 0.17677669529663687) # 1/(2**2.5) ~ 1/(2**1.5)
    ,'2nd_degree': (0.17677669529663687, 0.08838834764831843) # 1/(2**3.5) ~ 1/(2**2.5)
    ,'3rd_degree': (0.08838834764831843, 0.044194173824159216) # 1/(2**4.5) ~ 1/(2**3.5)
    ,'Unrelated':  (0.08838834764831843, 0) # < 1/(2**4.5)
}


# Original Code: Runpeng Luo
def kinship_inference_terry(ifile1: str, ifile2: str, ffile_i: str, jfile1: str, jfile2: str, ffile_j: str, arr_size: int):
    
    iid = decrypt_filename(ifile1)
    jid = decrypt_filename(jfile1)

    # both individuals are mapped to same reference, which leads to same SNP positions
    ifilter = get_pos_filter(ffile_i)
    jfilter = get_pos_filter(ffile_j)

    # only consider the SNP positions that neither individuals contain the missing alleles
    valid_filter = ifilter & jfilter
    # count common valid positions
    N_valid = valid_filter.count(1)
    print("Total collided valided SNP: ", N_valid)
    # obtain set intersections from raw data
    # inte_data_i/j contains the homozygote SNPs for i/j
    inte_data_i, homo_pos_i = process_raw(ifile1, ifile2, valid_filter)
    inte_data_j, homo_pos_j = process_raw(jfile1, jfile2, valid_filter)

    # homozygous positions bloom filter
    bloom_Hi = init_bloom_filter(arr_size, False, data_list=homo_pos_i, iid=iid)
    bloom_Hj = init_bloom_filter(arr_size, False, data_list=homo_pos_j, iid=jid)

    # set intersections from raw data bloom filter
    bloom_Sin = init_bloom_filter(arr_size, False, data_list=inte_data_i, iid=iid)
    bloom_Sjn = init_bloom_filter(arr_size, False, data_list=inte_data_j, iid=jid)

    # compute N^i_{Aa}
    Ni_Aa = N_valid - bloom_Sin.m
    Nj_Aa = N_valid - bloom_Sjn.m

    N_AaAa = N_valid - card_estim_union(bloom_Hi, bloom_Hj)
    N_AAaa = card_estim_intersection(bloom_Hi, bloom_Hj) - card_estim_intersection(bloom_Sin, bloom_Sjn)

    print("Ni_Aa: ", Ni_Aa)
    print("Nj_Aa: ", Nj_Aa)
    print("N_AaAa: ", N_AaAa)
    print("N_AAaa: ", N_AAaa)

    wtn_eps = within_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa)
    print("within family kinship eps: ", wtn_eps)
    get_relation(wtn_eps)

    btn_eps = between_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa)
    print("between family kinship eps: ", btn_eps)
    get_relation(btn_eps)

    return


def process_raw(file1: str, file2: str, pos_filter: bitarray):
    inte_data, homo_pos = [], []
    fd1, fd2 = open(file1, "r"), open(file2, "r")
    for valid_bit, line1, line2 in zip(pos_filter, fd1, fd2):
        if valid_bit != 1:
            continue
        pos, ref, alt1 = line1.strip().split(":")
        _,    _,  alt2 = line2.strip().split(":")

        if alt1 == alt2:
            inte_data.append(pos+":"+ref+":"+alt1)
            homo_pos.append(pos)
    return inte_data, homo_pos

def get_pos_filter(ffile: str) -> bitarray:
    fd = open(ffile, "r")
    pos_filter = bitarray(fd.readline()[:-1])
    fd.close()
    return pos_filter

def within_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa):
    numerator = N_AaAa - 2 * N_AAaa
    denominator = Ni_Aa + Nj_Aa
    return float(numerator) / denominator

def between_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa):
    fst_term = float(N_AaAa - 2 * N_AAaa) / (2 * Ni_Aa)
    snd_term = 0.5
    trd_term = (Ni_Aa + Nj_Aa) / float(4 * Ni_Aa)
    return fst_term + snd_term - trd_term

def get_relation(kin_eps: float):
    for relation, (up, low) in relationship_dict.items():
        if kin_eps <= up and kin_eps > low:
            print(relation)
            break
    

def gtruth(ifile1: str, ifile2: str, ffile_i: str, jfile1: str, jfile2: str, ffile_j: str):
    print(">>> Ground-truth solution: ")
    Ni_Aa = 0
    Nj_Aa = 0
    N_AaAa = 0
    N_AAaa = 0

    ifilter = get_pos_filter(ffile_i)
    jfilter = get_pos_filter(ffile_j)

    # only consider the SNP positions that neither individuals contain the missing alleles
    valid_filter = ifilter & jfilter
    valid_count = valid_filter.count(1)

    ifd1, ifd2, jfd1, jfd2 = [open(f, "r") for f in [ifile1, ifile2, jfile1, jfile2]]
    for vbit, il1, il2, jl1, jl2 in zip(valid_filter, ifd1, ifd2, jfd1, jfd2):
        if vbit != 1:
            continue
        i_homo = il1 == il2
        j_homo = jl1 == jl2
        if not i_homo:
            Ni_Aa += 1
        if not j_homo:
            Nj_Aa += 1
        if not i_homo and not j_homo:
            N_AaAa += 1
        if i_homo and j_homo and il1 != jl1:
            N_AAaa += 1
    ifd1.close()
    ifd2.close()
    jfd1.close()
    jfd2.close()

    print("N_valid: ", valid_count)
    print("Ni_Aa: ", Ni_Aa)
    print("Nj_Aa: ", Nj_Aa)
    print("N_AaAa: ", N_AaAa)
    print("N_AAaa: ", N_AAaa)

    wtn_eps = within_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa)
    print("within family kinship eps: ", wtn_eps)
    get_relation(wtn_eps)

    btn_eps = between_family_kinship_eps(Ni_Aa, Nj_Aa, N_AAaa, N_AaAa)
    print("between family kinship eps: ", btn_eps)
    get_relation(btn_eps)

    return