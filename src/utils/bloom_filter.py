import mmh3
from bitarray import bitarray

import math
import os

from utils.utils import *

# Original Code: Terry Yu
# Reformat and Optimize: John Luo

class Bloom_Filter:
    def __init__(self, iid, insertion_size: int, arr_size: int, k_size = 0, arr_str: str="") -> None:
        self.iid = iid
        self.m = insertion_size
        self.n = arr_size
        self.k = self.opt_hashcount() if k_size == 0 else k_size
        self.fp_prob1 = self.approx_fp_prob_v1()
        self.fp_prob2 = self.approx_fp_prob_v2()

        if arr_str != "":
            # from file
            self.bit_arr = bitarray(arr_str)
        else:
            # from scratch
            self.bit_arr = bitarray(self.n)
            self.bit_arr.setall(0)
        
        print(">>Bloom filter created, m: {0}, n: {1}, k: {2}, id: {3}".format(self.m, self.n, self.k, self.iid))
        print("FP rate (approx): {0}, FP rate (non-approx): {1}".format(self.fp_prob1, self.fp_prob2))
        print("Cardinality estimation: ", self.card_estim())
        return

    def opt_hashcount(self) -> int:
        return max(round(math.log(2) * float(self.n)/self.m), 1)
    
    def approx_fp_prob_v1(self) -> float:
        return math.pow(1 - math.exp(-self.k * self.m / float(self.n)), self.k)
    
    def approx_fp_prob_v2(self) -> float:
        return math.pow(1 - math.pow(1 - 1.0/self.n, self.k * self.m), self.k)

    def card_estim(self) -> int:
        num_ones = float(self.bit_arr.count(1))
        numerator = math.log(1 - num_ones / self.n)
        denominator = self.k * math.log(1 - 1.0 / self.n)
        return max(round(numerator / denominator), 0)

    def insert_element(self, elem) -> None:
        digests = []
        for i in range(self.k):
            digest = mmh3.hash(elem, i) % self.n
            digests.append(digest)
            self.bit_arr[digest] = 1
        return
    
    def query_element(self, elem) -> bool:
        for i in range(self.k):
            digest = mmh3.hash(elem, i) % self.n
            if self.bit_arr[digest] == 0:
                return False
        return True

    def get_intersect(self, bloom, set=False) -> bitarray:
        res = self.bit_arr & bloom.bit_arr
        if set:
            self.bit_arr = res
        return res
    
    def get_union(self, bloom, set=False):
        res = self.bit_arr | bloom.bit_arr
        if set:
            self.bit_arr = res
        return res

def init_bloom_filter(arr_size: int, from_file: bool, data_list=[], iid=None, filename=None):
    bloom = None
    if from_file:
        iid_f, isize = decrypt_filename(filename)
        bloom = Bloom_Filter(iid_f, isize, arr_size)
        print("Init from file: ", filename)
        with open(filename, "r") as fd:
            for line in fd:
                if line == "\n":
                    break
                bloom.insert_element(str(line[:-1]))
            fd.close()
    else:
        bloom = Bloom_Filter(iid, len(data_list), arr_size)
        print("Init from data list")
        for elem in data_list:
            bloom.insert_element(elem)
    print("Bloom filter initialized")
    return bloom

def store_bloom_tofile(bloom: Bloom_Filter, filename: str):
    os.system("echo "" > " + filename)
    with open(filename, "w") as fd:
        fd.write("{0}\t{1}\n".format(bloom.m, bloom.n))
        fd.write("{0}\n".format(bloom.bit_arr.to01()))
        fd.close()
    print("Store Bloom filter into file: ", filename)
    return 0

def get_bloom_fromfile(filename: str):
    fd = open(filename, "r")
    m, n = fd.readline()[:-1].split("\t")
    arr_str = fd.readline()[:-1]
    print("Obtained Bloom filter from file: ", filename)
    return Bloom_Filter(-1, int(m), int(n), str(arr_str))

    """Estimate cardinality between two bloom filters intersection
    """
def card_estim_intersection(bloom1: Bloom_Filter, bloom2: Bloom_Filter):
    card1 = bloom1.card_estim()
    card2 = bloom2.card_estim()
    temp_arr1 = bloom1.bit_arr
    bloom1.get_union(bloom2, True)
    card_union = bloom1.card_estim()
    bloom1.bit_arr = temp_arr1
    return card1 + card2 - card_union

    """Estimate cardinality between two bloom filters union
    """
def card_estim_union(bloom1: Bloom_Filter, bloom2: Bloom_Filter):
    temp_arr1 = bloom1.bit_arr
    bloom1.get_union(bloom2, True)
    card_union = bloom1.card_estim()
    bloom1.bit_arr = temp_arr1
    return card_union

def card_estim_triple_intersection():
    return

