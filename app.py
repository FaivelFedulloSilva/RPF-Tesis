from dash import Dash, dcc, html, Input, Output
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

import pandas as pd
import polars as pl
from polars import DataFrame

from file_handlers import GTFhandler, GTFobject, FastaHandler, BAMhandler
from libs import DNAStringSet

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
    

def detect_pauses(ts: list[int],min_distance, min_len):
    l = np.array(ts)
    [q1,q3] = np.quantile(l, [.25,.75])
    inter_quantile = q3-q1
    max = q3+inter_quantile
    return get_regions_from_outliers([[i, l[i]] for i in range(len(l)) if l[i] > max ],min_distance, min_len)

def get_pause_to_graph(pauses):
    traces = []
    for p in pauses:
        trace = [[],[]]
        trace[0].append(p['start'])        
        trace[1].append(p['values'][0])
        if p['start'] != p['end']:
            trace[0].append(p['end'])
            trace[1].append(p['values'][-1])
        traces.append(trace)
    return traces

#TODO: cambiar methods a una clase de contantes
def get_normalize_rpf_over_rna(cov_rpf: list[int], cov_rna: list[int], offset: int=1, window_size:int = 0, method: str= 'Offset'):
    if method == "Offset":
        new_rna = [c if c > 0 else offset for c in cov_rna]
    else:
        new_rna = [int(np.mean(cov_rna[i-window_size: i+window_size])) if cov_rna[i] == 0 else cov_rna[i]for i in range(window_size, len(cov_rna)-window_size)]
        new_rna = [c for c in cov_rna[:window_size]] + new_rna + [c for c in cov_rna[-window_size:]]
        new_rna = [c if c > 0 else offset for c in new_rna]
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

df_rpf: DataFrame = pl.read_json('coverage_rpf.json')
df_rna: DataFrame = pl.read_json('coverage_rna.json')

transcript_ids = df_rpf.select('transcript_id').to_series().to_list()





# dna = FastaHandler(REF_PATH)
# bam_rpf = BAMhandler(RPF_PATH, 'b')
# bam_rna = BAMhandler(TOTAL_PATH, 'b')


# gtf_file_handler = GTFhandler(GTF_PATH)
# gtf = gtf_file_handler.get_gtf_object()
# filtered_gtf = gtf.filter_by_feature('exon')
# filtered_gtf_transcripts = filtered_gtf.get_transcripts_ids()



app = Dash(__name__)

app.layout = html.Div([
    dcc.Graph(id='comparative-figure'),
    dcc.Dropdown(transcript_ids, value=transcript_ids[0], id='feature-select', ),
    dcc.Slider(0, 80, 1,
               value=10,
               id='my-slider'
    ),
    dcc.RadioItems(['Offset', 'Mean'], 'Offset', id='method'),
    dcc.Slider(1, 10, 1,
               value=3,
               id='window-size'
    ),
    html.H4("Maximum distance between two outliers to be considerer part of the same pause region "),
    dcc.Input(
        id="max-distance", type="number",
        min=2, max=30, step=1, value=2
    ),
    html.H4("Minimum number of outliers in row to be consider a pause region"),
    dcc.Input(
        id="min-length", type="number",
        min=1, max=100, step=1, value=3
    ),

    dcc.Graph(id="graph"),
])


@app.callback(
    Output('comparative-figure', 'figure'),
    Output('graph', 'figure'),
    Input('feature-select', 'value'),
    Input('window-size', 'value'),
    Input('method', 'value'),
    Input('max-distance', 'value'),
    Input('min-length', 'value'),\
)
def update_figure(selected_feature, window_size, method, max_distance, min_len):
    
    # df = filtered_gtf.get_transcripts_data([selected_feature]) 

    # feature_sequences_rpf = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rpf)
    # feature_sequences_rna = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rna)

    feature_sequences_rpf = df_rpf.filter(pl.col('transcript_id') == selected_feature)
    feature_sequences_rna = df_rna.filter(pl.col('transcript_id') == selected_feature)

    # cov_rpf = feature_sequences_rpf[selected_feature][1]  
    # cov_rna = feature_sequences_rna[selected_feature][1]

    cov_rpf = feature_sequences_rpf.select('coverage_per_base').to_series().to_list()[0]
    cov_rna = feature_sequences_rna.select('coverage_per_base').to_series().to_list()[0]
    cov_rpf_over_rna = get_normalize_rpf_over_rna(cov_rpf, cov_rna, offset=1,window_size=window_size,method=method)

    outliers = detect_pauses(cov_rpf_over_rna, max_distance, min_len)
    pauses = get_pause_to_graph(outliers)

    print(pauses)

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
        go.Bar(x=position, y=cov_rpf_over_rna), row=3, col=1
        )
    for pause in pauses:
        fig.add_trace(
            go.Scatter(x=pause[0], y=pause[1],line = dict(color = ('rgb(160, 160, 160)')),showlegend=False), row=3, col=1,
        )


    fig2 = make_subplots(rows=1, cols=2, shared_xaxes=True,subplot_titles=(f"RPF boxplot" f"RPF/RNA boxplot"))
    fig2.add_trace(go.Box(y=cov_rpf, quartilemethod="linear"), row=1, col=1)
    fig2.add_trace(go.Box(y=cov_rpf_over_rna, quartilemethod="linear"), row=1, col=2)

    return fig,fig2

@app.callback(
    Output('feature-select', 'options'),
    Input('my-slider', 'value')
)
def update_dropdown(selected_filter):
    transcript_ids_filtered = (df_rpf
                .lazy()
                .filter(pl.col('coverage_percentage') > selected_filter/100)
                .select('transcript_id')
            ).collect().to_series().to_list()

    return transcript_ids_filtered



if __name__ == '__main__':
    app.run_server(debug=True)
