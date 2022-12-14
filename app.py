from dash import Dash, dcc, html, Input, Output
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

import pandas as pd
from polars import DataFrame

from file_handlers import GTFhandler, GTFobject, FastaHandler, BAMhandler
from libs import DNAStringSet

def get_normalize_rpf_over_rna(cov_rpf: list[int], cov_rna: list[int]):
    new_rna = [c if c > 0 else 1 for c in cov_rna]
    return [cov_rpf[i]/new_rna[i] for i in range(len(cov_rpf))]

def get_second_elements(l: list[str,str])->list[int]:
    new_list = []
    for elem in l:
        new_list.append(int(elem[1]))

    return new_list

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


RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'
REF_PATH = 'Data/reference/hg38.fa'
GTF_PATH = 'Data/genesFiltrada.gtf'

dna = FastaHandler(REF_PATH)
bam_rpf = BAMhandler(RPF_PATH, 'b')
bam_rna = BAMhandler(TOTAL_PATH, 'b')


gtf_file_handler = GTFhandler(GTF_PATH)
gtf = gtf_file_handler.get_gtf_object()
filtered_gtf = gtf.filter_by_feature('exon')
filtered_gtf_transcripts = filtered_gtf.get_transcripts_ids()

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminderDataFiveYear.csv')

app = Dash(__name__)

app.layout = html.Div([
    dcc.Graph(id='comparative-figure'),
    dcc.Dropdown(filtered_gtf_transcripts, value=filtered_gtf_transcripts[0], id='feature-select', ),
])


@app.callback(
    Output('comparative-figure', 'figure'),
    Input('feature-select', 'value')
)
def update_figure(selected_feature):
    df = filtered_gtf.get_transcripts_data([selected_feature]) 

    feature_sequences_rpf = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rpf)
    feature_sequences_rna = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rna)

    cov_rpf = feature_sequences_rpf[selected_feature][1]  
    cov_rna = feature_sequences_rna[selected_feature][1]

    position = np.arange(0, len(cov_rpf))
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,subplot_titles=(f"RPF for {selected_feature}", f"RNA for {selected_feature}", f"RPF/RNA for {selected_feature}"))
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

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
