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

def get_base_sequence(gtf: DataFrame, dna: DNAStringSet = None):
    print(gtf.sample(20))
    q = (gtf
            .lazy()
            .groupby(by='transcript_id')
            .agg(
                [pl.col('seqname').unique().first(), pl.col('strand'), pl.col('start'),  pl.col('end')]
                )
        )
    newdf = q.collect().sample(5)
    print(newdf)
    iter_df = newdf.iterrows(named=True)
    feature_sequences = {}
    for row in iter_df:
        feature_seq = []
        ref_seq = dna.get_sequence(row.seqname)
        for index in range(len(row.start)):
            seq = ref_seq.get_sub_sequence(row.start[index], row.end[index])
            if row.strand[index] == '-':
                seq = seq.reverse_complement()
            feature_seq.append(seq.get_as_string())
        feature_sequences[row.transcript_id] = ''.join(feature_seq)
    for key, value in feature_sequences.items():
        print(key, " -> ", value)

gtf = GTFhandler(GTF_PATH)
filtered_gtf = gtf.filter_by_feature('exon')
dna = FastaHandler(REF_PATH)
get_base_sequence(filtered_gtf, dna.get_reference())