import numpy as np
import pysam
from file_handlers import GTFhandler
from file_handlers import FastaHandler
from file_handlers import BAMhandler
from file_handlers import GTFobject
import pandas as pd
import polars as pl
from polars import DataFrame
from libs.DNAStringSet import DNAStringSet
from collections import deque
from math import trunc

import plotly.graph_objects as go
import plotly.express as px


from timeit import default_timer as timer

from libs.DNAString import DNAString

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'
REF_PATH = 'Data/reference/hg38.fa'
GTF_PATH = 'Data/genesFiltrada.gtf'

def get_regions_from_outliers(l: list[list[int,int]], min_distance: int, min_len: int):
    regions = [[l[0]]]
    last_outlier = l[0]
    for outlier in l[1:]:
        if outlier[0] - last_outlier[0] < min_distance:
            regions[-1].append(outlier)
        else:
            regions.append([outlier])
        last_outlier = outlier
    filtered_regions = [region for region in regions if len(region) >= min_len ]

    parsed_regions =[]
    region = {'start': float('inf'), 'end': float('-inf'), 'values': []}
    for reg in filtered_regions:
        for r in reg:
            if r[0] < region['start']:
                region['start'] = r[0]
            if r[0] > region['end']:
                region['end'] = r[0]
            region['values'].append(r[1])
        parsed_regions.append(region)
        region = {'start': float('inf'), 'end': float('-inf'), 'values': []}

    return parsed_regions
    

def detect_outliers(ts: list[int],min_distance, min_len):
    l = np.array(ts)
    [q1,q3] = np.quantile(l, [.25,.75])
    inter_quantile = q3-q1
    max = q3+inter_quantile
    print(max)
    return get_regions_from_outliers([[i, l[i]] for i in range(len(l)) if l[i] > max ],min_distance, min_len)

# TODO Actualizar para que se peuda pasar como parametro la forma de calcular. Por ejemplo usando una media movil con minimo de 1 para rna
def get_normalize_rpf_over_rna(cov_rpf: list[int], cov_rna: list[int]):
    new_rna = [c if c > 0 else 1 for c in cov_rna]
    return [cov_rpf[i]/new_rna[i] for i in range(len(cov_rpf))]

def plot_comparative(cov_rpf: list[int], cov_rna: list[int], name: str):
    from plotly.subplots import make_subplots
    position = np.arange(0, len(cov_rna))
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,subplot_titles=(f"RPF for {name}", f"RNA for {name}", f"RPF/RNA for {name}"))
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


def filter_coverage(cov, filter = 0):
    """ Return true if the quantity of non-zero values in cov is more than 100-filter percent"""
    # return np.count_nonzero(cov)/np.size(cov) >= (100-filter)/100
    return np.count_nonzero(cov)/np.size(cov)


def get_feature_sequence_and_coverage(df: DataFrame, dna: DNAStringSet, bam: BAMhandler, filter: int):
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
        # if filter_coverage(np.array(feature_depth), filter):
        feature_sequences[row.transcript_id] = [''.join(feature_seq), filter_coverage(feature_depth) , feature_depth]
    return feature_sequences

def get_feature_sequence_and_coverage_and_save_to_file(df: DataFrame, dna: DNAStringSet, bam: BAMhandler, file_path: str):
    iter_df = df.iterrows(named=True)
    feature_sequences = {
        'transcript_id': [],
        'coverage_percentage': [],
        'sequence': [],
        'coverage_per_base': []
        }
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
        # if filter_coverage(np.array(feature_depth), filter):
        # feature_sequences[row.transcript_id] = [''.join(feature_seq), filter_coverage(feature_depth), feature_depth]
        feature_sequences['transcript_id'].append(row.transcript_id)
        feature_sequences['coverage_percentage'].append(filter_coverage(feature_depth))
        feature_sequences['sequence'].append("".join(feature_seq))
        feature_sequences['coverage_per_base'].append(feature_depth)
    dt = pl.DataFrame(feature_sequences)
    dt.write_json(file_path)

# df = pl.read_json('coverage_rpf.json')
# print(df)
# 
# input()
gtf_file_handler = GTFhandler(GTF_PATH)
gtf = gtf_file_handler.get_gtf_object()
filtered_gtf = gtf.filter_by_feature('exon')
dna = FastaHandler(REF_PATH, polars=False)
bam_rpf = BAMhandler(RPF_PATH, 'b')
bam_rna = BAMhandler(TOTAL_PATH, 'b')

# print(gtf.get_transcripts_ids())
transcripts = filtered_gtf.get_transcripts_ids()
# print(transcripts)
# print(transcripts[0].item())

# input()

for x in [10, 100, 1000]:
# newdf = get_base_sequence(filtered_gtf, dna.get_reference(), bam_rpf, bam_rna)
    df = filtered_gtf.get_transcripts_data(transcripts[0:x]) 
# feature_sequences_rpf = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rpf, 80)

    start = timer()
    get_feature_sequence_and_coverage_and_save_to_file(df, dna.get_reference(), bam_rpf, f'coverage_rpf_{x}.json')

    end = timer()
    print(end - start)
# print(feature_sequences_rpf)
# feature_sequences_rna = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rna)
# print(rollling_function(feature_sequences_rpf['NM_033657'][1],50,0,0))
# plot_comparative(feature_sequences_rpf['NM_001127322'][1],feature_sequences_rna['NM_001127322'][1], 'NM_001127322')
# plot_comparative(feature_sequences_rpf['NM_033657'][1],feature_sequences_rna['NM_033657'][1], 'NM_033657')