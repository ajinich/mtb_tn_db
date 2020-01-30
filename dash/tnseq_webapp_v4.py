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
import seaborn as sns
import visdcc
from numpy import inf

external_stylesheets = [dbc.themes.UNITED]
path_data = '../data/'
main_data = pd.read_csv('Tn_library_DASH.csv', header=0)
# main_data = main_data.drop(columns='gene_name')
column_descriptors = pd.read_csv('column_descriptors.csv', header=0)
genes_data = pd.merge(main_data, column_descriptors,
                      left_on='Expt', right_on='column_ID', how='left')
genes_data = genes_data.drop(
    columns=['sequence_data', 'column_ID', 'journal', 'first_author', 'last_author'])

cogs_df = pd.read_csv('all_cogs.csv', header=0)
cogs_desc = pd.read_csv('cog_names.csv', header=0, index_col=0, squeeze=True)
orf_details = pd.read_csv('ORF_details_final.csv')
unique_expts = main_data['Expt'].unique()
unique_Rvs = list(main_data['Rv_ID'].unique())
unique_genes = list(main_data['gene_name'].unique())
main_data['id'] = main_data['Rv_ID']
main_data.set_index('id', inplace=True, drop=False)

df_uk = pd.read_csv(os.path.join(
    path_data, 'unknown_essentials/unknown_ALL_levels_essential_scores.csv'))
df_uk = df_uk[['Rv_ID', 'gene_name', 'UK_score_4']]


def Table(dataframe, id):
    paper_URLs = dataframe['paper_URL']
    dataframe = dataframe.drop(columns=['paper_URL'])
    rows = []
    for i in range(len(dataframe)):
        sig = False
        row = []
        for col in dataframe.columns:
            val = dataframe.iloc[i][col]
            # update this depending on which
            # columns you want to show links for
            # and what you want those links to be
            if col == 'q-val':
                if val <= 0.05:
                    sig = True
            if col == 'paper_title':
                link = paper_URLs.iloc[i]
                cell = html.Td(
                    html.A(href=link, target="_blank", children=val))
            else:
                cell = html.Td(children=val)
            row.append(cell)
        if sig:
            rows.append(html.Tr(row, style={'background-color': '#eeeeee'}))
        else:
            rows.append(html.Tr(row))
    Thead = html.Thead(
        [html.Tr([html.Th(col, style={'color': "#DD4814", 'text-align': 'center',
                                      'background-color': '#cccccc'}) for col in dataframe.columns])]
    )
    Tbody = html.Tbody(rows, style={'text-align': 'center'})
    # columns=html.Col()
    return html.Table([Thead, Tbody], id=id, className="display", style={'border': '1', 'border-collapse': 'separate'})


def discretize_q_values(row):
    q_val = row['q-val']
    if q_val < 0.01:
        q_val_d = 3
    elif q_val < 0.05:
        q_val_d = 2
    else:
        q_val_d = 1
    return q_val_d


def unknown_essential_xy(TnSeq_screen, df_data, df_uk, rand_param=0.6):
    # Grab data for a single TnSeq screen
    df_data_test = main_data[main_data['Expt'] == TnSeq_screen]
    df_data_test['q_val_D'] = df_data_test.apply(discretize_q_values, 1)

    # Merge with unknowns:
    df_vis = df_data_test.merge(df_uk, on=['Rv_ID', 'gene_name'], how='inner')

    # Get x-y datasets:
    rv_ids = df_vis.Rv_ID.values
    uk_list = np.array(df_vis.UK_score_4)
    q_list = np.array(df_vis.q_val_D)

    # randomize:
    uk_rd = np.array([uk + rand_param * random() -
                      rand_param / 2 for uk in uk_list])
    q_rd = np.array([q + rand_param * random() -
                     rand_param / 2 for q in q_list])

    # color the unknown-essentials differently:
    current_palette = sns.color_palette()
    # all genes are gray by default.
    color_list = np.array([(0.35, 0.35, 0.35)] * df_vis.shape[0])
    # Unknown essentials in a different color.
    ind_temp = list(df_vis[(df_vis.q_val_D == 3) &
                           (df_vis.UK_score_4 == 4)].index)
    color_list[ind_temp] = current_palette[0]
    color_list_rgb = ['rgb(' + ', '.join([str(np.round(rgb, 2))
                                          for rgb in col]) + ')' for col in color_list]

    return uk_rd, q_rd, color_list_rgb, rv_ids


external_scripts = ['https://code.jquery.com/jquery-3.3.1.min.js',
                    'https://cdn.datatables.net/v/dt/dt-1.10.18/datatables.min.js']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                external_scripts=external_scripts)

navbar = dbc.NavbarSimple([
    dbc.NavItem(dbc.NavLink('Analyze datasets',
                            active=True, href='/analyze_datasets')),
    dbc.NavItem(dbc.NavLink('Analyze genes', href='/analyze_genes')),
    dbc.NavItem(dbc.NavLink('About', active=True, href='/about'))
], brand="Mtb Tn-seq database", color='primary', light=True)

analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id='Sel_dataset',
                                                  options=[{'label': x, 'value': x}
                                                           for x in unique_expts],
                                                  value=unique_expts[0])
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
                                 style={'background-color': '#f5f5f5', 'padding': '30px', 'border-radius': '25px',
                                        'border-color': '#dcdcdc', 'border-width': '2px', 'border-style': 'solid'}),
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
                                 ], align='center', width=5),
                                 dbc.Col([
                                     html.Div([
                                         html.Label('Gene List')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=3),
                             ]),
                             dbc.Row([
                                 dbc.Col([
                                     # html.Div([
                                     #     html.Label('About this dataset')
                                     # ], style={'textAlign': 'center', 'display': 'block'}),
                                     html.Div(id='dataset_metadata')], width=3,
                                     style={'background-color': '#f5f5f5', 'padding': '30px', 'border-radius': '25px',
                                            'border-color': '#dcdcdc', 'border-width': '2px', 'border-style': 'solid'}),
                                 dbc.Col([
                                     dcc.Graph(id='volcano'),
                                 ],
                                     width=5, align='center'),
                                 dbc.Col([
                                     dt.DataTable(id='table',
                                                  columns=[{"name": i, "id": i} for i in [
                                                      'Rv_ID', 'log2FC', 'q-val']],
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
                                                  style_as_list_view=True,
                                                  tooltip={
                                                      "Rv_ID": "this is a test tooltip"}
                                                  )
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
                                     dcc.Dropdown(id='Sel_cog',
                                                  options=[{'label': x, 'value': x} for x in
                                                           ['Under-represented', 'Over-represented']],
                                                  value='Under-represented'),
                                     dcc.Graph(id='cog')

                                 ], width=6, align='center'),
                                 dbc.Col([
                                     dcc.Graph(id='tSNE_plot')
                                 ], width=6, align='center')
                             ], justify='center')
                             ])

analyze_genes = html.Div([
    dbc.Row([html.Label('Pick a gene')]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='Sel_gene', options=[{'label': x, 'value': x} for x in unique_genes+unique_Rvs],
                         placeholder='Select a gene', multi=False)]),
        dbc.Col([
            html.Div(id='orf_metadata')])
    ], style={'background-color': '#f5f5f5', 'padding': '30px', 'border-radius': '25px',
              'border-color': '#dcdcdc', 'border-width': '2px', 'border-style': 'solid'}),
    html.Br(),
    html.Br(),
    html.Div(id='dataset_table_div')
])

about = html.Div([
    html.Label('Desc and credits'),
    html.Br(),
    html.A('Download raw data',
           href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # dbc.Button("Download raw data", href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # html.Label('Download raw_data')

])

# app.layout = html.Div([navbar, body])
app.layout = html.Div(
    [
        dcc.Location(id="url", pathname="/analyze_datasets", refresh=False),
        navbar,
        dbc.Container(id="content", style={"padding": "20px"}),
        # You need to write this here for dash to load visdcc across the whole app
        visdcc.Run_js(id='javascript', run="$('#datatable1').DataTable()")
    ])

app.config.suppress_callback_exceptions = True
app.scripts.config.serve_locally = True


# TODO: still get javascript error with Rv1784, possible input df issue
@app.callback(dash.dependencies.Output("dataset_table_div", "children"),
              [dash.dependencies.Input("Sel_gene", 'value')]
              )
def generate_dataset_table(selected_gene):
    if selected_gene in unique_Rvs:
        df = genes_data[genes_data['Rv_ID'] == selected_gene]
    elif selected_gene in unique_genes:
        df = genes_data[genes_data['gene_name'] == selected_gene]
    else:
        return None
    df['q-val'] = np.round(df['q-val'], 2)
    df['log2FC'] = np.round(df['log2FC'], 2)
    df = df.sort_values(by='q-val')
    df = df[['Rv_ID', 'gene_name', 'paper_URL', 'paper_title', 'log2FC',
             'q-val', 'in_vitro_media', 'carbon_source', 'stat_analysis']]
    return [visdcc.Run_js(id='javascript', run="$('#datatable1').DataTable()"), Table(df, id='datatable1')]


@app.callback(dash.dependencies.Output("content", "children"), [dash.dependencies.Input("url", "pathname")])
def display_page(pathname):
    if pathname == "/analyze_datasets":
        return analyze_datasets
    if pathname == "/analyze_genes":
        return analyze_genes
    if pathname == '/about':
        return about
    # if not recognised, return 404 message
    return html.P("404 - page not found")


@app.callback([
    dash.dependencies.Output('download_dataset', 'href'),
    dash.dependencies.Output('download_dataset', 'download')
],
    [dash.dependencies.Input('Sel_dataset', 'value')])
def update_download_dataset(sel_dataset):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    csv_string = selected_data.to_csv(encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + \
        urllib.parse.quote(csv_string)
    download_string = sel_dataset + '.csv'
    return csv_string, download_string


@app.callback(
    dash.dependencies.Output('volcano', 'figure'),
    [dash.dependencies.Input('Sel_dataset', 'value'),
     dash.dependencies.Input('log2FC', 'value'),
     dash.dependencies.Input('q-val', 'value'),
     dash.dependencies.Input('table', "derived_virtual_row_ids"),
     dash.dependencies.Input('table', "selected_row_ids")])
def update_figure(sel_dataset, log2FC, pval, row_ids, selected_row_ids):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    if selected_row_ids is None:
        selected_row_ids = []

    if row_ids is None:
        dff = selected_data
        row_ids = selected_data['id']
    else:
        dff = selected_data.loc[row_ids]
    # print('here', dff['q-val'])
    # print ('here', np.unique(-np.log10(dff['q-val'])))
    max_log_pval = np.unique(-np.log10(dff['q-val']))[-2]
    # print(max_log_pval)
    inf_repl = np.ceil(max_log_pval) + 1
    # print(inf_repl)
    dff['pval_plotting'] = -np.log10(dff['q-val'])
    dff['pval_plotting'].replace(np.inf, inf_repl, inplace=True)
    # print('here2', np.arange(0, inf_repl + 0.5, 0.5))
    tickvals = list(np.arange(0, inf_repl + 0.5, 0.5))
    ticklab = tickvals.copy()
    ticklab[-1] = 'Inf'

    for_x_ticks = dff['log2FC']
    for_x_ticks = for_x_ticks.replace([np.inf, -np.inf], np.nan)
    for_x_ticks = for_x_ticks.dropna()
    tickvals_x = list(np.arange(int(for_x_ticks.min() - 1),
                                int(for_x_ticks.max() + 1), 1))
    ticklab_x = tickvals_x.copy()

    ticked = dff.index.isin(selected_row_ids)
    ticked_data = dff[ticked]
    unticked_data = dff[~ticked]
    generated_filter = (unticked_data['q-val'] <= pval) & (
        (unticked_data['log2FC'] <= (-log2FC)) | (unticked_data['log2FC'] >= log2FC))
    sig_data = unticked_data[generated_filter]
    non_sig_data = unticked_data[~generated_filter]

    traces = []
    traces.append(go.Scatter(x=sig_data['log2FC'],
                             y=sig_data['pval_plotting'],
                             text=sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Outside cutoff',
                             marker={'opacity': 0.6, 'size': 10,
                                     'color': 'orangered'},
                             showlegend=False,
                             ))
    traces.append(go.Scatter(x=non_sig_data['log2FC'],
                             y=non_sig_data['pval_plotting'],
                             text=non_sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Pass cutoff',
                             marker={'opacity': 0.6,
                                     'size': 10, 'color': 'grey'},
                             showlegend=False
                             ))
    traces.append(go.Scatter(x=ticked_data['log2FC'],
                             y=ticked_data['pval_plotting'],
                             text=ticked_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers+text',
                             textposition='bottom center',
                             name='T',
                             marker={'opacity': 0.6,
                                     'size': 10, 'color': 'green'},
                             showlegend=False
                             ))
    return {'data': traces,
            'layout': go.Layout(
                autosize=False,
                margin={
                    'l': 45,
                    'r': 15,
                    'pad': 0,
                    't': 30,
                    'b': 90, },
                xaxis={'title': 'log2FC', 'ticktext': ticklab_x,
                       'tickvals': tickvals_x},
                yaxis={
                    'title': '-log10(q-val)', 'ticktext': ticklab, 'tickvals': tickvals},
                hovermode='closest'
            )}


@app.callback(
    dash.dependencies.Output('table', 'data'),
    [dash.dependencies.Input('Sel_dataset', 'value')])
def update_datatable(sel_dataset):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    selected_data = selected_data.drop('Expt', axis=1)
    selected_data['q-val'] = np.round(selected_data['q-val'], 2)
    selected_data['log2FC'] = np.round(selected_data['log2FC'], 2)
    selected_data = selected_data.sort_values(by='log2FC')
    return selected_data.to_dict('records')


@app.callback(
    dash.dependencies.Output('tSNE_plot', 'figure'),
    [dash.dependencies.Input('Sel_dataset', 'value')])
def update_bubble(sel_dataset):
    uk_rd, q_rd, color_list_rgb, rv_ids = unknown_essential_xy(
        sel_dataset, main_data, df_uk)
    return {
        'data': [
            go.Scatter(
                x=uk_rd,
                y=q_rd,
                text=rv_ids,
                mode='markers',
                opacity=1,
                hoverinfo='text',
                marker={
                    'size': 15,
                    'line': {'width': 1.5, 'color': 'black'},
                    'color': color_list_rgb
                }
            )
        ],
        'layout': go.Layout(
            autosize=True,
            # width=800,
            # height=500,
            xaxis=go.layout.XAxis(
                tickmode='array',
                tickvals=[0, 1, 2, 3, 4],
                ticktext=['most well<br>characterized', '', '', '',
                          'least well<br>characterized'],
                tickfont=dict(size=14),
                title='Annotation'

            ),
            yaxis=go.layout.YAxis(
                tickmode='array',
                tickvals=[1, 2, 3],
                ticktext=['non-essential', 'q-val < 0.05', 'q-val < 0.01'],
                tickangle=270,
                tickfont=dict(size=14),
                title='Essentiality'
            ),
            margin={'l': 30, 'b': 100, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest'
        )
    }


@app.callback(dash.dependencies.Output('cog', 'figure'),
              [dash.dependencies.Input('Sel_dataset', 'value'),
               dash.dependencies.Input('Sel_cog', 'value'),
               dash.dependencies.Input('log2FC', 'value'),
               dash.dependencies.Input('q-val', 'value')])
def update_cog(sel_dataset, sel_cog, log2FC, pval):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    if sel_cog == 'Under-represented':
        sel_subset_filter = (
            selected_data['q-val'] <= pval) & (selected_data['log2FC'] <= -log2FC)
        colorscale = 'Cividis'
    else:
        sel_subset_filter = (
            selected_data['q-val'] <= pval) & (selected_data['log2FC'] >= -log2FC)
        colorscale = 'Viridis'
    sel_subset = selected_data[sel_subset_filter]
    cog_total_freq = cogs_df['COG'].value_counts(normalize=True)
    sel_cogs = cogs_df[cogs_df['X.Orf'].isin(sel_subset['Rv_ID'])]
    sel_cogs_freq = sel_cogs['COG'].value_counts(normalize=True)
    normalized_cogs = sel_cogs_freq / cog_total_freq
    normalized_cogs = normalized_cogs[~normalized_cogs.isnull()]
    normalized_cogs = normalized_cogs.sort_values()
    cog_names = cogs_desc.loc[list(normalized_cogs.index)]
    cog_names = list(cog_names.values)

    bar_data = [go.Bar(y=list(normalized_cogs.index), x=list(normalized_cogs.values),
                       orientation='h',
                       text=cog_names,
                       hoverinfo='text',
                       marker={'color': list(normalized_cogs.values), 'colorscale': colorscale})]
    return {'data': bar_data,
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
            }


@app.callback(
    dash.dependencies.Output(component_id='orf_metadata',
                             component_property='children'),
    [dash.dependencies.Input(component_id='Sel_gene', component_property='value')])
def print_orf_metadata(sel_gene):
    sel_details = orf_details[orf_details['ORF'] == sel_gene]
    return sel_details['Description']


if __name__ == '__main__':
    app.run_server(debug=True)