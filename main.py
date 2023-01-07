import numpy as np
import matplotlib as plt
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


import numpy as np

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, RangeTool
from bokeh.plotting import figure, show, output_file, save

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'
REF_PATH = 'Data/reference/hg38.fa'
GTF_PATH = 'Data/genesFiltrada.gtf'


def plot_with_plotly(cov: list[int]):
    print(cov)
    position = np.arange(0, len(cov))
    covdf = {
        'location': position,
        'coverage': cov   
    }
    covdf = pd.DataFrame(covdf)

    fig = px.bar(covdf, x='location', y='coverage', title='Time Series with Rangeslider')

    fig.update_xaxes(rangeslider_visible=True)
    fig.show()




def get_second_elements(l: list[str,str])->list[int]:
    new_list = []
    for elem in l:
        new_list.append(int(elem[1]))

    return new_list

def get_base_sequence(gtf: DataFrame, dna: DNAStringSet = None, bam: BAMhandler = None):
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
        feature_depth = []
        ref_seq = dna.get_sequence(row.seqname)
        for index in range(len(row.start)):
            seq = ref_seq.get_sub_sequence(row.start[index], row.end[index])
            if row.strand[index] == '-':
                seq = seq.reverse_complement()
            print(row.seqname, row.start[index], row.end[index])
            cov = bam.region_coverage(row.seqname, row.start[index], row.end[index])
            feature_depth += get_second_elements(cov)
            feature_seq.append(seq.get_as_string())
        feature_sequences[row.transcript_id] = [''.join(feature_seq), feature_depth]
    for key, value in feature_sequences.items():
        plot_with_plotly(value[1])
        input()
        print(key, " -> ", value)

gtf = GTFhandler(GTF_PATH)
filtered_gtf = gtf.filter_by_feature('exon')
dna = FastaHandler(REF_PATH)
bam = BAMhandler(RPF_PATH, 'b')
get_base_sequence(filtered_gtf, dna.get_reference(), bam)