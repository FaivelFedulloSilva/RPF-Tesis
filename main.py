import numpy as np
import pysam
from file_handlers import GTFhandler
from file_handlers import FastaHandler
from file_handlers import BAMhandler
import pandas as pd
import polars as pl
from polars import DataFrame
from libs.DNAStringSet import DNAStringSet

import plotly.graph_objects as go
import plotly.express as px

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'
REF_PATH = 'Data/reference/hg38.fa'
GTF_PATH = 'Data/genesFiltrada.gtf'


def get_normalize_rpf_over_rna(cov_rpf: list[int], cov_rna: list[int]):
    new_rna = [c if c > 0 else 1 for c in cov_rna]
    return [cov_rpf[i]/new_rna[i] for i in range(len(cov_rpf))]

def plot_comparative(cov_rpf: list[int], cov_rna: list[int], name: str):
    from plotly.subplots import make_subplots
    position = np.arange(0, len(cov_rna))
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True)
    fig.update_yaxes(fixedrange=True)
    fig.add_trace(
        go.Bar(x=position, y=cov_rpf), row=1, col=1
        )   
    fig.add_trace(
        go.Bar(x=position, y=cov_rna), row=2, col=1
        )
    fig.add_trace(
        go.Bar(x=position, y=get_normalize_rpf_over_rna(cov_rpf, cov_rna)), row=3, col=1
        )
    fig.show()

def plot_with_plotly(cov: list[int], name: str):
    print(cov)
    position = np.arange(0, len(cov))
    covdf = {
        'location': position,
        'coverage': cov   
    }
    covdf = pd.DataFrame(covdf)

    fig = px.bar(covdf, x='location', y='coverage', title=f'Coverage per base of transcript {name}')

    fig.update_xaxes(rangeslider_visible=True)
    fig.show()




def get_second_elements(l: list[str,str])->list[int]:
    new_list = []
    for elem in l:
        new_list.append(int(elem[1]))

    return new_list


def get_lazy_query(transcript: str, gtf: DataFrame):
    query = (gtf
            .lazy()
            .groupby(by='transcript_id')
            .agg(
                [pl.col('seqname').unique().first(), pl.col('strand'), pl.col('start'),  pl.col('end')]
                )
            .filter(pl.col('transcript_id') == transcript)
        )

    return query

def get_feature_sequence_and_coverage(df: DataFrame, dna: DNAStringSet, bam: BAMhandler):
    iter_df = df.iterrows(named=True)
    feature_sequences = {}
    for row in iter_df:
        feature_seq = []
        feature_depth = []
        ref_seq = dna.get_sequence(row.seqname)
        for index in range(len(row.start)):
            seq = ref_seq.get_sub_sequence(row.start[index], row.end[index])
            if row.strand[index] == '-':
                seq = seq.reverse_complement()
            # print(row.seqname, row.start[index], row.end[index])
            cov = bam.region_coverage(row.seqname, row.start[index], row.end[index])
            feature_depth += get_second_elements(cov)
            feature_seq.append(seq.get_as_string())
        feature_sequences[row.transcript_id] = [''.join(feature_seq), feature_depth]
    return feature_sequences

def get_base_sequence(gtf: DataFrame, dna: DNAStringSet = None, bam_rpf: BAMhandler = None, bam_rna: BAMhandler = None):
    # print(gtf.sample(20))
    # q = (gtf
    #         .lazy()
    #         .groupby(by='transcript_id')
    #         .agg(
    #             [pl.col('seqname').unique().first(), pl.col('strand'), pl.col('start'),  pl.col('end')]
    #             )
    #         .filter(pl.col('transcript_id') == 'NM_001127322')
    #     )
    q = get_lazy_query('NM_001127322', gtf)
    r = get_lazy_query('NM_033657', gtf)
    newdf = q.collect()
    newdf.extend(r.collect())
    print(newdf)
    feature_sequences_rpf = get_feature_sequence_and_coverage(newdf, dna, bam_rpf)
    feature_sequences_rna = get_feature_sequence_and_coverage(newdf, dna, bam_rna)
    plot_comparative(feature_sequences_rpf['NM_033657'][1],feature_sequences_rna['NM_033657'][1], 'NM_033657')
    input()
    for key, value in feature_sequences_rpf.items():
        plot_with_plotly(value[1], key)
        input()
        print(key, " -> ", value)
    for key, value in feature_sequences_rna.items():
        plot_with_plotly(value[1], key)
        input()
        print(key, " -> ", value)

gtf = GTFhandler(GTF_PATH)
filtered_gtf = gtf.filter_by_feature('exon')
dna = FastaHandler(REF_PATH)
bam_rpf = BAMhandler(RPF_PATH, 'b')
bam_rna = BAMhandler(TOTAL_PATH, 'b')
get_base_sequence(filtered_gtf, dna.get_reference(), bam_rpf, bam_rna)