import numpy as np
import matplotlib as plt
import pysam
from file_handlers import GTFhandler
from file_handlers import FastaHandler
import polars as pl
from polars import DataFrame
from libs.DNAStringSet import DNAStringSet

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'
REF_PATH = 'Data/reference/hg38.fa'
GTF_PATH = 'Data/genesFiltrada.gtf'

def get_base_sequence(gtf: DataFrame, dna: DNAStringSet):
    pass


gtf = GTFhandler(GTF_PATH)
filtered_gtf = gtf.filter_by_feature('exon')
dna = FastaHandler(REF_PATH)
get_base_sequence(filtered_gtf, dna.get_reference())