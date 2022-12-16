import numpy as np
import matplotlib as plt
import pysam

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'

rpf_bam = pysam.AlignmentFile(RPF_PATH, 'rb')

iter = rpf_bam.fetch("chr1")
iter = list(iter)
print(iter[10].get_forward_sequence())
print(iter[10].get_reference_sequence())
