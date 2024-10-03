import os
import urllib.parse
from io import StringIO
from random import random

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import dash_table as dt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from numpy import inf
import itertools

# TODO: FIX BUBBLE PLOT FROM NOTEBOOK AND CHANGE XAXIS LABELS!!!

#####
# SECTION 1: Read in data, create static variables
#####
external_stylesheets = [dbc.themes.UNITED]
path_data = '../data/'

#read in data
std_data = pd.read_csv(os.path.join(
    path_data, 'standardized_data_dash.tsv'), sep='\t', dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': np.float, 'q-val': np.float})
si_data = pd.read_csv(os.path.join(
    path_data, 'si_data_dash.tsv'), sep='\t', dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': np.float, 'q-val': np.float})
gene_metadata_df = pd.read_csv(os.path.join(
    path_data, 'gene_metadata_dash.tsv'), sep='\t')
# TODO: make num_replicates into int
# TODO: fill empty strings in meaning etc with ' ' instead of nan
metadata = pd.read_csv(os.path.join(path_data, 'col_desc_dash.tsv'), sep='\t')

# make dictionaries that lookup si_data from std_data and vice versa
dict_std_to_si = dict(zip(metadata.column_ID_std, metadata.column_ID_SI))
dict_si_to_std = dict(zip(metadata.column_ID_SI, metadata.column_ID_std))
# make dictionary that tells you if SI data has enough datapoints to be plotted
dict_plot_si = dict(zip(metadata.column_ID_SI, metadata.plot_SI_graph))

# make lists of unique expts/rvids/genes for dropdown lists
unique_expts = list(metadata['column_ID_std'].unique()) + \
    list(metadata['column_ID_SI'].unique())
unique_expts = [x for x in unique_expts if str(x) != 'nan']
unique_Rvs = sorted(gene_metadata_df.Rv_ID.to_list())
unique_genes = sorted(gene_metadata_df.gene_name.to_list())
unique_genes = [x for x in unique_genes if x != '-']

# do data wrangling on si and std data for dash requirements
std_data['id'] = std_data['Rv_ID']
std_data.set_index('id', inplace=True, drop=False)
si_data['id'] = si_data['Rv_ID']
si_data.set_index('id', inplace=True, drop=False)

# read in cog annotations, cogs_desc should be read as a pd.Series
path_annotation = '../data/annotations/'
cogs_df = pd.read_csv(os.path.join(path_annotation, 'all_cogs.csv'))
cogs_desc = pd.read_csv(os.path.join(
    path_annotation, 'cog_names.csv'), header=0, index_col=0, squeeze=True)


# select columns of dataset metadata to show in analyse genes table
metadata_cols_display = [col for col in metadata.columns if col not in [
    'plot_SI_graph', 'column_ID_std', 'column_ID_SI', 'paper_title', 'paper_URL']]
sel_gene_table_columns = [{"name": i,
                           "id": i,
                           } for i in [
    'Rv_ID', 'gene_name', 'Expt', 'log2FC', 'q-val'] + metadata_cols_display]

# show paper column as markdown so that html links work
sel_gene_table_columns.append(
    {"name": 'paper', "id": 'paper', "presentation": 'markdown'})

# list plotly buttons to remove from all graphs
plotly_buttons_remove = [
    'pan2d', 'lasso2d', 'select2d', 'hoverClosestCartesian', 'hoverCompareCartesian',
    'zoomIn2d', 'zoomOut2d', 'toggleSpikelines', 'autoScale2d', 'zoom2d'
]


def empty_plot(label_annotation):
    """Returns a dictionary with elements for an empty plot with centered text

    Args:
        label_annotation (str): Text to display

    Returns:
        dict: Elements for plotting
    """
    return {
        "layout": {
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "annotations": [
                {"text": label_annotation,
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {"size": 16}
                 }]
        }}


def unknown_essential_xy(selected_data):
    """Generate plotting data for bubble plot

    Args:
        selected_data (pd.DataFrame): Data filtered on selected experiment

    Returns:
        Multiple outputs required for plotting
    """

    df_vis = selected_data.dropna(subset=['annotation_score'])
    df_vis = df_vis.reset_index(drop=True)

    # NO MORE RANDOMIZATION! Instead this:
    rv_id_list = []
    x_coords_list = []
    y_coords_list = []
    color_list = []
    scatter_size_list = []
    lw_list = []

    for ann in [1, 2, 3, 4, 5]:
        for qq in [1, 2, 3]:

            df = df_vis[(df_vis.q_val_D.values == qq) &
                        (df_vis.annotation_score.values == ann)]
            # update rv_id list
            rv_id_list.append(list(df.Rv_ID.values))

            if df.shape[0] < 30:
                scatter_size = 6
                edge_param = 0.40
                lw = 4
            elif df.shape[0] < 100:
                scatter_size = 6
                edge_param = 0.40
                lw = 3
            elif df.shape[0] < 200:
                scatter_size = 3
                edge_param = 0.45
                lw = 2
            elif df.shape[0] < 400:
                scatter_size = 3
                edge_param = 0.45
                lw = 2
            else:
                scatter_size = 3
                edge_param = 0.45
                lw = 1

            # Update scatter marker size
            scatter_size_list += [scatter_size] * df.shape[0]
            # Update line-width size
            lw_list += [lw] * df.shape[0]

            if ann <= 2 and qq >= 2:
                color_temp = '#2b7bba'
            else:
                color_temp = '#585858'

            # update color_list
            color_list += [color_temp] * df.shape[0]

            # num_sqrs = int(np.ceil(np.sqrt(df.shape[0])))
            num_sqrs = np.max([int(np.ceil(np.sqrt(df.shape[0]))), 10])

            xrange = np.linspace(ann - edge_param, ann + edge_param, num_sqrs)
            yrange = np.linspace(qq + edge_param, qq - edge_param, num_sqrs)

            coords = list(itertools.product(yrange, xrange))
            coords = coords[:df.shape[0]]
            x_coords = [c[1] for c in coords]
            y_coords = [c[0] for c in coords]
            # update coordinates
            x_coords_list += x_coords
            y_coords_list += y_coords

    return x_coords_list, y_coords_list, color_list, rv_id_list, scatter_size_list, lw_list


def filter_dataset(sel_dataset, sel_standardized):
    """Using outputs from sel_dataset and sel_standardized return appropriate data and name

    Args:
        sel_dataset (str): Selected dataset eg: griffin_glycerol_vs_mbio_H37Rv
        sel_standardized (str): Selected standardized eg: 'Original'

    Returns:
        pd.DataFrame: Filtered data
        str: Dataset name.
    """
    if sel_standardized == 'Standardized':
        dataset_name = dict_si_to_std.get(sel_dataset, sel_dataset)
        dff = std_data[std_data['Expt'] == dataset_name]
    else:
        dataset_name = dict_std_to_si.get(sel_dataset, sel_dataset)
        dff = si_data[si_data['Expt'] == dataset_name]
    return dff, dataset_name


def update_download_dataset(dff, dataset_name):
    """Given data filtered by experiment, provide data for downloading

    Args:
        dff (pd.DataFrame): Data filtered by experiment
        dataset_name (str): Dataset name

    Returns:
        str: download string
        str: download file name
    """
    dff = dff.copy(deep=True)
    dff.reset_index(inplace=True, drop=True)
    dff = dff.drop(columns='id')
    csv_string = dff.to_csv(encoding='utf-8', sep='\t', index=False)
    csv_string = "data:text/plain;charset=utf-8," + \
        urllib.parse.quote(csv_string)
    download_string = dataset_name + '.tsv'
    return csv_string, download_string


# TODO: decide what to do with NAs
def update_num_significant(dff, log2FC, qval):
    """Given filtered data, log2FC and qval return text indicating num significant genes

    Args:
        dff (pd.DataFrame): Data filtered by expt
        log2FC (float): log2FC cutoff
        qval (float): qval cutoff

    Returns:
        list: List of html elements as text 
    """
    dff = dff.dropna(axis=0, subset=['log2FC', 'q-val'])
    num_neg_sig = dff[(dff['q-val'] <= qval) &
                      (dff['log2FC'] <= -log2FC)].shape[0]
    num_pos_sig = dff[(dff['q-val'] <= qval) &
                      (dff['log2FC'] >= log2FC)].shape[0]
    text = [html.Span(f'Number of neg_sig : {num_neg_sig} ; Number of pos_sig : {num_pos_sig}'),
            ]
    return text


#####
# SECTION 2: LAYOUT
#####
# initialize app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# create navbar
navbar = dbc.NavbarSimple([
    dbc.NavItem(dbc.NavLink('Analyze datasets',
                            active=True, href=app.get_relative_path('/analyze_datasets'))),
    dbc.NavItem(dbc.NavLink('Analyze genes',
                            href=app.get_relative_path('/analyze_genes'))),
    dbc.NavItem(dbc.NavLink('About', active=True,
                            href=app.get_relative_path('/about')))
], brand="Mtb Tn-seq database", color='primary', light=True)

# Layout for page analyze datasets
analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id='sel_dataset',
                                                  options=[{'label': x, 'value': x}
                                                           for x in unique_expts],
                                                  value=unique_expts[0], clearable=False),
                                     dcc.Dropdown(
                                         id='sel_standardized', clearable=False)
                                 ], width=4),
                                 dbc.Col([
                                     daq.Slider(id='log2FC', min=0, max=6, value=1, step=0.5,
                                                size=300,
                                                marks={x: x for x in range(0, 7)}, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "log2FC"})
                                 ], width=4),
                                 dbc.Col([
                                     daq.Slider(id='q-val', min=0, max=1, value=0.05, step=0.05,
                                                marks={
                                                    x / 10: x / 10 for x in range(1, 11)},
                                                size=300, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "q-val"})
                                 ], width=4),
                                 html.Br(),
                                 html.A('Download this dataset', id='download_dataset', download="", href="",
                                        target="_blank"),

                             ], align='center',
    style={'background-color': '#f5f5f5',
                                 'padding': '30px',
                                 'border-radius': '25px',
                                 'border-color': '#dcdcdc',
                                 'border-width': '2px',
                                 'border-style': 'solid'}),
    html.Br(),
    html.Br(),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Label('About this dataset')
            ], style={'textAlign': 'center', 'display': 'block'}),
        ], align='center', width=3),
        dbc.Col([
            html.Div([
                html.Label('Volcano plot')
            ], style={'textAlign': 'center', 'display': 'block'}),
            html.Div(id='num_significant', style={
                     'textAlign': 'center', 'display': 'block'}),
        ], align='center', width=5),
        dbc.Col([
            html.Div([
                html.Label('Gene List')
            ], style={'textAlign': 'center', 'display': 'block'}),
        ], align='center', width=3),
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='dataset_metadata')], width=3,
            style={'background-color': '#f5f5f5',
                   'padding': '30px',
                   'border-radius': '25px',
                   'border-color': '#dcdcdc',
                   'border-width': '2px',
                   'border-style': 'solid'}),
        dbc.Col([dcc.Loading(id='loading_volcano', children=dcc.Graph(id='volcano')),
                 ],
                width=5, align='center'),
        dbc.Col([dcc.Loading(id='loading_dataset_table', children=dt.DataTable(id='sel_dataset_table',
                                                                               columns=[{"name": i, "id": i} for i in [
                                                                                   'Rv_ID', 'gene_name', 'log2FC', 'q-val']],
                                                                               sort_action='native',
                                                                               row_selectable='multi',
                                                                               selected_rows=[],
                                                                               page_action='native',
                                                                               page_size=15,
                                                                               page_current=0,
                                                                               style_header={'color': '#e95420', 'font-weight': 'bold',
                                                                                             'text-align': 'center'},
                                                                               style_cell_conditional=[
                                                                                   {'if': {'column_id': 'q-val'},
                                                                                    'width': '30%'}
                                                                               ],
                                                                               style_data_conditional=[
                                                                                   {'if': {'row_index': 'odd'},
                                                                                    'backgroundColor': 'rgb(248,248,248)'}
                                                                               ],
                                                                               style_cell={
                                                                                   'font-family': 'ubuntu',
                                                                                   'font-size': 14,
                                                                                   'height': '10px',
                                                                                   'textOverFlow': 'ellipsis',
                                                                                   'text-align': 'center',
                                                                                   'overflow': 'hidden'
                                                                               },
                                                                               #    style_as_list_view=True,
                                                                               ))
                 ], width=4, align='center')
    ]),
    html.Br(),
    html.Br(),
    html.Br(),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Label('COG Categories')
            ], style={'textAlign': 'center', 'display': 'block'}),
        ], align='center', width=6),
        dbc.Col([
            html.Div([
                html.Label('Essentiality plot')
            ], style={'textAlign': 'center', 'display': 'block'}),
        ], align='center', width=6),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='sel_cog',
                         options=[{'label': x, 'value': x} for x in
                                  ['Under-represented', 'Over-represented']],
                         value='Under-represented'),
            dcc.Graph(id='cog')

        ], width=6, align='center'),
        dbc.Col([dcc.Loading(id='loading_bubble', children=dcc.Graph(id='bubble_plot'))
                 ], width=6, align='center')
    ], justify='center')
])

# Layout for page analyze genes
analyze_genes = html.Div([
    dbc.Row([html.Label('Pick a gene')]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='sel_gene',
                         options=[{'label': x, 'value': x}
                                  for x in unique_genes + unique_Rvs],
                         placeholder='Select a gene',
                         multi=False,
                         searchable=True),
            dcc.Dropdown(id='sel_standardized_gene_table',
                         options=[
                             {'label': x, 'value': x} for x in ['Standardized', 'Original']],
                         value='Standardized', clearable=False,
                         multi=False)
        ]),
        dbc.Col([
            html.Div(id='gene_metadata')])
    ], style={'background-color': '#f5f5f5',
              'padding': '30px',
              'border-radius': '25px',
              'border-color': '#dcdcdc',
              'border-width': '2px',
              'border-style': 'solid'}),
    html.Br(),
    html.Br(),
    dt.DataTable(id='sel_gene_table',
                 columns=sel_gene_table_columns,
                 sort_action='native',
                 page_action='native',
                 #  filter_action='native',
                 page_size=15,
                 page_current=0,
                 style_header={'color': '#e95420', 'font-weight': 'bold',
                               'textAlign': 'center'},
                 style_data_conditional=[
                     {'if': {'row_index': 'odd'},
                      'backgroundColor': 'rgb(248,248,248)'}
                 ],
                 style_cell={
                     'font-family': 'ubuntu',
                     'font-size': 14,
                     'height': '10px',
                     # 'textOverFlow': 'ellipsis',
                     'textAlign': 'center',
                     'overflow': 'hidden'
                 })

])

# Layout for page About
about = html.Div([
    html.H5('Contact'),
    html.Span("For bug reports and data submissions, contact "),
    html.A('Adrian Jinich',
           href="mailto:adj2010@med.cornell.edu", target='_blank'),
    html.Br(),
    html.Br(),
    html.H5('Raw data'),
    html.Span('Raw data is available '),
    html.A('here',
           href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # dbc.Button("Download raw data", href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # html.Label('Download raw_data')
])

# app.layout dynamically takes in different content based on the path. See next callback
app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        navbar,
        dbc.Container(id="content", style={"padding": "20px"}),
    ])

app.config.suppress_callback_exceptions = True
app.scripts.config.serve_locally = True

#####
# SECTION 3: CALLBACKS
#####


@ app.callback(Output("content", "children"),
               [Input("url", "pathname")])
def display_content(path):
    """
    Takes in path from the URL and returns layout for one of three pages
    """
    page_name = app.strip_relative_path(path)
    if page_name == 'analyze_datasets':
        return analyze_datasets
    if page_name == "analyze_genes":
        return analyze_genes
    if page_name == 'about':
        return about


@ app.callback(
    [Output('sel_standardized', 'options'),
     Output('sel_standardized', 'value')],
    [Input('sel_dataset', 'value')])
def update_standardized_dropdown(sel_dataset):
    """Take in the dataset selected and update the sel_standardized dropdown

    Args:
        sel_dataset (str): Selected dataset eg: griffin_glycerol_vs_mbio_H37Rv

    Returns:
        list: List of options for sel_standardized
        str: Default value to display from the options
    """
    # is it a std_dataset?
    if sel_dataset in dict_std_to_si:
        # does a corresponding original (aka si) dataset exist?
        if pd.isna(dict_std_to_si[sel_dataset]):
            return [{'label': x, 'value': x} for x in ['Standardized']], 'Standardized'
        else:
            return [{'label': x, 'value': x} for x in ['Standardized', 'Original']], 'Standardized'
    # is it an original dataset?
    else:
        # does a corresponding std_dataset exist?
        if pd.isna(dict_si_to_std[sel_dataset]):
            return [{'label': x, 'value': x} for x in ['Original']], 'Original'
        else:
            return [{'label': x, 'value': x} for x in ['Standardized', 'Original']], 'Original'


@ app.callback([
    Output('download_dataset', 'href'),
    Output('download_dataset', 'download'),
    Output('num_significant', 'children'),
],
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value'),
     Input('log2FC', 'value'),
     Input('q-val', 'value'),
     ])
def update_multiple_outputs_analyze_datasets(sel_dataset, sel_standardized, log2FC, qval):
    """Using user inputs, update download dataset and num significant

    Args:
        sel_dataset (str): User selected dataset 
        sel_standardized (str): User selected standardized/original
        log2FC (float): User selected log2FC cutoff
        qval (flaot): User selected qval cutoff

    Returns:
        str: href for download
        str: file name for download
        list: list of text for number of significant. 
              ' ' is returned if not enough genes in experiment for this to be meaningful
    """
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    num_significant_text = update_num_significant(dff, log2FC, qval)
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            num_significant_text = ' '
    csv_string, download_string = update_download_dataset(dff, dataset_name)

    return csv_string, download_string, num_significant_text


@app.callback([Output('volcano', 'figure'),
               Output('volcano', 'config')],
              [Input('sel_dataset', 'value'),
               Input('sel_standardized', 'value'),
               Input('log2FC', 'value'),
               Input('q-val', 'value'),
               Input('sel_dataset_table', "derived_virtual_selected_row_ids")
               ])
def update_volcano(sel_dataset, sel_standardized, log2FC, qval, selected_row_ids):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {
        'modeBarButtonsToRemove': plotly_buttons_remove,
        'toImageButtonOptions': {
            'height': 700,
            'width': 700,
            'scale': 5,
            'filename': f'{dataset_name}_volcano.png'
        }
    }
    if sel_standardized == 'Original':
        # Is there enough data for a meaningful plot?
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    # weird plotly requirement
    if selected_row_ids is None:
        selected_row_ids = []
    # TODO: Figure out NAs
    dff = dff.dropna(axis=0, subset=['log2FC', 'q-val'])
    # make qval ticks, replacing the np.nans with inf
    # what is current second highest max log10 transformed qval?
    # note that max will be inf
    max_log_qval = np.unique(-np.log10(dff['q-val']))[-2]
    # create new column in dff with qval for plotting, replace inf values
    inf_repl = np.ceil(max_log_qval) + 1
    dff['qval_plotting'] = -np.log10(dff['q-val'])
    dff['qval_plotting'].replace(np.inf, inf_repl, inplace=True)

    # create x and y tick vals and labels
    tickvals = list(np.arange(0, inf_repl + 0.5, 0.5))
    ticklab = tickvals.copy()
    ticklab[-1] = 'Inf'
    for_x_ticks = dff['log2FC']
    # TODO: WHAT is this? - commented for now
    # for_x_ticks = for_x_ticks.replace([np.inf, -np.inf], np.nan)
    # for_x_ticks = for_x_ticks.dropna()
    tickvals_x = list(np.arange(int(for_x_ticks.min() - 1),
                                int(for_x_ticks.max() + 1), 1))
    ticklab_x = tickvals_x.copy()

    # split data into selected (ie ticked), unselected_sig, unselected_non_sig
    ticked = dff['id'].isin(selected_row_ids)
    ticked_data = dff[ticked]
    unticked_data = dff[~ticked]
    generated_filter = (unticked_data['q-val'] <= qval) & (
        (unticked_data['log2FC'] <= (-log2FC)) | (unticked_data['log2FC'] >= log2FC))
    sig_data = unticked_data[generated_filter]
    non_sig_data = unticked_data[~generated_filter]

    # make traces for each kind of data
    traces = []
    traces.append(go.Scatter(x=sig_data['log2FC'],
                             y=sig_data['qval_plotting'],
                             text=sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Outside cutoff',
                             marker={'opacity': 0.6, 'size': 10,
                                       'color': 'orangered'},
                             showlegend=False,
                             ))
    traces.append(go.Scatter(x=non_sig_data['log2FC'],
                             y=non_sig_data['qval_plotting'],
                             text=non_sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Pass cutoff',
                             marker={'opacity': 0.6,
                                       'size': 10,
                                       'color': 'grey'},
                             showlegend=False
                             ))
    traces.append(go.Scatter(x=ticked_data['log2FC'],
                             y=ticked_data['qval_plotting'],
                             text=ticked_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers+text',
                             textposition='bottom center',
                             name='T',
                             marker={'opacity': 0.6,
                                       'size': 10,
                                       'color': 'green'},
                             showlegend=False
                             ))
    # return dict of plotting and config for plotly
    return ({'data': traces,
             'layout': go.Layout(
                 autosize=False,
                 margin={'l': 45, 'r': 15, 'pad': 0, 't': 30, 'b': 90},
                 xaxis={'title': 'log2FC', 'ticktext': ticklab_x,
                        'tickvals': tickvals_x, 'fixedrange': True},
                 yaxis={'title': '-log10(q-val)', 'ticktext': ticklab,
                        'tickvals': tickvals, 'fixedrange': True},
                 hovermode='closest'
             )}, config)


@ app.callback(
    [Output('bubble_plot', 'figure'),
     Output('bubble_plot', 'config')],
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value'),
     ])
def update_bubble(sel_dataset, sel_standardized):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {'modeBarButtonsToRemove': plotly_buttons_remove, 'toImageButtonOptions': {
        'height': 500,
        'width': 700, 'scale': 5, 'filename': f'{dataset_name}_bubble.png'
    }}
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    x_coords_list, y_coords_list, color_list, rv_id_list, scatter_size_list, lw_list = unknown_essential_xy(
        dff)
    print(scatter_size_list)
    return ({
        'data': [
            go.Scatter(
                x=x_coords_list,
                y=y_coords_list,
                text=rv_id_list,
                mode='markers',
                # opacity=1,
                hoverinfo='text',
                marker_size=scatter_size_list,
                marker={
                    # 'size': 15,
                    'line': {'width': 1.5, 'color': 'black'},
                    'color': color_list
                }
            )
        ],
        'layout': go.Layout(
            autosize=True,
            # width=800,
            # height=500,
            xaxis=go.layout.XAxis(
                tickmode='array',
                tickvals=[1, 2, 3, 4, 5],
                ticktext=['most well<br>characterized', '', '', '',
                          'least well<br>characterized'],
                tickfont=dict(size=14),
                title='Annotation',
                showgrid=False

            ),
            yaxis=go.layout.YAxis(
                tickmode='array',
                tickvals=[1, 2, 3],
                ticktext=['non-essential', 'q-val < 0.05', 'q-val < 0.01'],
                tickangle=270,
                tickfont=dict(size=14),
                title='Essentiality',
                showgrid=False
            ),
            margin={'l': 30, 'b': 100, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest'
        )
    }, config)


@ app.callback(
    Output('dataset_metadata', 'children'),
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value')])
def print_dataset_metadata(sel_dataset, sel_standardized):
    if sel_standardized == 'Standardized':
        dataset_name = dict_si_to_std.get(sel_dataset, sel_dataset)
        dff = metadata[metadata['column_ID_std'] == dataset_name]
    else:
        dataset_name = dict_std_to_si.get(sel_dataset, sel_dataset)
        dff = metadata[metadata['column_ID_SI'] == dataset_name]
    text = [
        html.Strong('Summary'),
        html.Span(': ' + dff['meaning'].values[0]),
        html.Br(),
        html.Br(),
        html.Strong('Original Publication:'),
        html.Br(),
        html.A(dff['paper_title'].values[0],
               href=dff['paper_URL'].values[0], target='_blank'),
        html.Br(),
        html.Br(),
        html.Strong('No of control replicates'),
        html.Span(': ' + str(dff['num_replicates_control'].values[0])),
        html.Br(),
        html.Br(),
        html.Strong('No of experimental replicates'),
        html.Span(': ' + str(dff['num_replicates_experimental'].values[0]))
    ]
    return text


@ app.callback(
    Output('sel_dataset_table', 'data'),
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value')])
def update_dataset_table(sel_dataset, sel_standardized):
    dff, _ = filter_dataset(sel_dataset, sel_standardized)
    dff = dff[['Rv_ID', 'gene_name', 'log2FC', 'q-val', 'id']]
    dff['q-val'] = dff['q-val'].astype(float).round(2)
    dff['log2FC'] = dff['log2FC'].astype(float).round(2)
    dff = dff.sort_values(by='log2FC')
    return dff.to_dict('records')


@ app.callback([Output('cog', 'figure'), Output('cog', 'config')],
               [Input('sel_dataset', 'value'),
                Input('sel_standardized', 'value'),
                Input('sel_cog', 'value'),
                Input('log2FC', 'value'),
                Input('q-val', 'value')])
def update_cog(sel_dataset, sel_standardized, sel_cog, log2FC, qval):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {'modeBarButtonsToRemove': plotly_buttons_remove, 'toImageButtonOptions': {
        'height': 500,
        'width': 700, 'scale': 5, 'filename': f'{dataset_name}_bubble.png'
    }}
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    if sel_cog == 'Under-represented':
        sel_subset_filter = (
            dff['q-val'] <= qval) & (dff['log2FC'] <= -log2FC)
        colorscale = 'Cividis'
    else:
        sel_subset_filter = (
            dff['q-val'] <= qval) & (dff['log2FC'] >= -log2FC)
        colorscale = 'Viridis'
    sel_subset = dff[sel_subset_filter]
    # calculate genome freq
    cog_total_freq = cogs_df['COG'].value_counts(normalize=True)
    # calculate subset freq
    sel_cogs = cogs_df[cogs_df['X.Orf'].isin(sel_subset['Rv_ID'])]
    sel_cogs_freq = sel_cogs['COG'].value_counts(normalize=True)
    # calculate enrichment
    normalized_cogs = sel_cogs_freq / cog_total_freq
    # format
    normalized_cogs = normalized_cogs[~normalized_cogs.isnull()]
    normalized_cogs = normalized_cogs.sort_values()
    cog_names = cogs_desc.loc[list(normalized_cogs.index)]
    cog_names = list(cog_names.values)

    bar_data = [go.Bar(y=list(normalized_cogs.index), x=list(normalized_cogs.values),
                       orientation='h',
                       text=cog_names,
                       hoverinfo='text',
                       marker={'color': list(normalized_cogs.values), 'colorscale': colorscale})]
    return ({'data': bar_data,
             'layout': go.Layout(
                 margin={
                     'l': 50,
                     'r': 10,
                     'pad': 3,
                     't': 30,
                     'b': 90, },
                 # paper_bgcolor='rgba(0,0,0,0)',
                 # plot_bgcolor = 'rgba(0,0,0,0)',
                 xaxis={'title': 'Normalized to genomic frequency'},
                 hovermode='closest',
                 shapes=[{'type': 'line', 'x0': 1, 'y0': 0, 'x1': 1, 'y1': len(normalized_cogs),
                          'line': {'color': 'grey', 'width': 1, 'dash': 'dot'}}])
             }, config)


@ app.callback(
    Output('sel_gene_table', 'data'),
    [Input('sel_gene', 'value'),
     Input('sel_standardized_gene_table', 'value')])
def update_genes_table(selected_gene, sel_standardized_gene_table):
    if sel_standardized_gene_table == 'Standardized':
        dff = std_data.copy()
        metadata_col = 'column_ID_std'
    else:
        dff = si_data.copy()
        metadata_col = 'column_ID_SI'
    if selected_gene in unique_Rvs:
        dff = dff[dff['Rv_ID'] == selected_gene]
    elif selected_gene in unique_genes:
        dff = dff[dff['gene_name'] == selected_gene]
    else:
        raise PreventUpdate
    metadata['paper'] = '[' + metadata['paper_title'] + \
        '](' + metadata['paper_URL'] + ')'
    metadata_cols_display = [col for col in metadata.columns if col not in [
        'plot_SI_graph', 'column_ID_std', 'column_ID_SI', 'paper_title', 'paper_URL']]
    metadata_trunc = metadata[[metadata_col] + metadata_cols_display]
    metadata_trunc = metadata_trunc.rename(columns={metadata_col: 'Expt'})
    merged_data = dff.merge(metadata_trunc, how='left', on='Expt')
    merged_data['q-val'] = np.round(merged_data['q-val'], 2)
    merged_data['log2FC'] = np.round(merged_data['log2FC'], 2)
    merged_data = merged_data.sort_values(by='q-val')
    return merged_data.to_dict('records')


@ app.callback(
    Output('gene_metadata', 'children'),
    [Input('sel_gene', 'value')])
def print_gene_metadata(sel_gene):
    if sel_gene in unique_Rvs:
        sel_details = gene_metadata_df[gene_metadata_df['Rv_ID'] == sel_gene]
    elif sel_gene in unique_genes:
        sel_details = gene_metadata_df[gene_metadata_df['gene_name'] == sel_gene]
        # sel_details = si_data[si_data['gene_name'] == sel_gene]
    else:
        return "gene not found"
    text = [
        html.Span(list(sel_details['Description'])[0]),
        html.Br(),
        html.Strong('mBio Call: '),
        html.Span(list(sel_details['Final Call'])[0])
    ]
    return text


if __name__ == '__main__':
    app.run_server(debug=True)
