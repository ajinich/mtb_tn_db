import dash
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_table as dt
import pandas as pd
import numpy as np
from numpy import inf
import plotly.graph_objs as go

external_stylesheets = [dbc.themes.UNITED]

main_data = pd.read_csv('Tn_library_DASH.csv', header=0)
main_data = main_data.drop(columns='gene_name')
cogs_df = pd.read_csv('all_cogs.csv', header=0)
cogs_desc = pd.read_csv('cog_names.csv', header=0, index_col=0, squeeze=True)
orf_details=pd.read_csv('ORF_details_final.csv')
unique_expts = main_data['Expt'].unique()
unique_genes = main_data['Rv_ID'].unique()


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

navbar = dbc.NavbarSimple([
    dbc.NavItem(dbc.NavLink('Analyze datasets', active=True, href='/analyze_datasets')),
    dbc.NavItem(dbc.NavLink('Analyze genes', href='/analyze_genes'))
], brand="Mtb Tn-seq database", color='primary', light=True)

analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id='Sel_dataset',
                                                  options=[{'label': x, 'value': x} for x in unique_expts],
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
                                                marks={x / 10: x / 10 for x in range(1, 11)},
                                                size=300, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "q-val"})
                                 ], width=4)
                             ], align='center',
                                 style={'background-color': '#f5f5f5', 'padding': '30px', 'border-radius': '25px',
                                        'border-color': '#dcdcdc', 'border-width': '2px', 'border-style': 'solid'}),
                             html.Br(),
                             html.Br(),
                             dbc.Row([
                                 dbc.Col([
                                     dt.DataTable(id='table',
                                                  columns=[{"name": i, "id": i} for i in ['Rv_ID', 'log2FC', 'q-val']],
                                                  n_fixed_rows=1,
                                                  sorting=True,
                                                  sorting_type='multi',
                                                  #sorting_settings=[],
                                                  row_selectable='multi',
                                                  style_table={
                                                      'maxHeight': 600,
                                                      # 'overflowY': 'auto'
                                                  },
                                                  style_header={'color': '#e95420', 'font-weight': 'bold',
                                                                'text-align': 'center'},
                                                  style_cell_conditional=[
                                                      {'if': {'column_id': 'q-val'},
                                                       'width': '30%'}
                                                      # {'if': {'row_index': 'odd'},
                                                      # 'backgroundColor': 'rgb(248,248,248)'}
                                                  ],
                                                  style_cell={
                                                      'font-family': 'ubuntu', 'font-size': 14, 'text-align': 'center'},
                                                  pagination_settings={'current_page': 0, 'page_size': 15,
                                                                    'displayed_pages': 1},
                                                  pagination_mode='fe',
                                                  navigation='page',
                                                  selected_rows=[],
                                                  style_as_list_view=True
                                                  )
                                 ], width=4),
                                 dbc.Col([
                                     html.Label('Volcano plot'),
                                     dcc.Graph(id='volcano'),
                                 ],
                                     width=4),
                                 dbc.Col([
                                     html.Label('COG categories'),
                                     dcc.Dropdown(id='Sel_cog',
                                                  options=[{'label': x, 'value': x} for x in
                                                           ['Under-represented', 'Over-represented']],
                                                  value='Under-represented'),
                                     dcc.Graph(id='cog')

                                 ], width=4)
                             ])])

analyze_genes = html.Div([
    dbc.Row([
        dbc.Col([
            html.Label('Pick a gene'),
            html.Br(),
            dcc.Dropdown(id='Sel_gene', options=[{'label': x, 'value': x} for x in unique_genes],
                         value='Rv0001', multi=False),
            html.Br(),
            html.Div(id='orf_metadata')
            # html.Br(),
            # html.Label('test'),
            # html.Br(),
            # html.Label('test2')
        ], width=4, style={'background-color': '#f5f5f5', 'padding': '30px', 'border-radius': '25px',
                           'border-color': '#dcdcdc', 'border-width': '2px', 'border-style': 'solid'}),
        dbc.Col([
            dt.DataTable(id='dataset_table',
                         columns=[{"name": i, "id": i} for i in main_data.columns],
                         sorting=True,
                         style_header={'color': '#e95420', 'font-weight': 'bold'},
                         style_cell_conditional=[
                             # {'if': {'row_index': 'odd'},
                             # 'backgroundColor': 'rgb(248,248,248)'},
                             {'if': {'column_id': 'Expt'},
                              'width': '40%'},
                             {'if': {'column_id': 'Rv_ID'},
                              'width': '20%'}
                         ],
                         style_cell={
                             'font-family': 'ubuntu', 'font-size': 14, 'text-align': 'center'},
                         pagination_settings={'current_page': 0, 'page_size': 15, 'displayed_pages': 1},
                         pagination_mode='fe',
                         navigation='page',
                         # selected_rows=[],
                         # style_as_list_view=True
                         )
        ], width=8)
    ])
])

# app.layout = html.Div([navbar, body])
app.layout = html.Div(
    [
        dcc.Location(id="url", pathname="/analyze_datasets"),
        navbar,
        dbc.Container(id="content", style={"padding": "20px"})
    ])

app.config.suppress_callback_exceptions = True


@app.callback(dash.dependencies.Output("content", "children"), [dash.dependencies.Input("url", "pathname")])
def display_page(pathname):
    if pathname == "/analyze_datasets":
        return analyze_datasets
    if pathname == "/analyze_genes":
        return analyze_genes
    # if not recognised, return 404 message
    return html.P("404 - page not found")


@app.callback(
    dash.dependencies.Output('volcano', 'figure'),
    [dash.dependencies.Input('Sel_dataset', 'value'),
     dash.dependencies.Input('log2FC', 'value'),
     dash.dependencies.Input('q-val', 'value'),
     dash.dependencies.Input('table', "derived_virtual_data"),
     dash.dependencies.Input('table', "derived_virtual_selected_rows")])
def update_figure(sel_dataset, log2FC, pval, rows, derived_virtual_selected_rows):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    if derived_virtual_selected_rows is None:
        derived_virtual_selected_rows = []

    if rows is None:
        dff = selected_data
    else:
        dff = pd.DataFrame(rows)
    # print('here', dff['q-val'])
    max_log_pval = np.unique(-np.log10(dff['q-val']))[-2]
    # print(max_log_pval)
    inf_repl = np.ceil(max_log_pval) + 1
    # print(inf_repl)
    dff['pval_plotting'] = -np.log10(dff['q-val'])
    dff['pval_plotting'].replace(inf, inf_repl, inplace=True)
    tickvals = list(np.arange(0, inf_repl + 0.5, 0.5))
    ticklab = tickvals.copy()
    ticklab[-1] = 'Inf'

    ticked = dff.index.isin(derived_virtual_selected_rows)
    ticked_data = dff[ticked]
    unticked_data = dff[~ticked]
    generated_filter = (unticked_data['q-val'] <= pval) & (
            (unticked_data['log2FC'] <= (-log2FC)) | (unticked_data['log2FC'] >= log2FC))
    sig_data = unticked_data[generated_filter]
    non_sig_data = unticked_data[~generated_filter]

    #
    #    colors=[]
    #    for i in range(len(dff)):
    #        if i in derived_virtual_selected_rows:
    #            colors.append('green')
    #        elif generated_filter[i]:
    #            colors.append('orange')
    #        else:
    #            colors.append('grey')
    #
    #    if label_selected:
    #        mode_for_selected='markers+text'
    #    else:
    #        mode_for_selected='markers'

    traces = []
    traces.append(go.Scatter(x=sig_data['log2FC'],
                             y=sig_data['pval_plotting'],
                             text=sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Outside cutoff',
                             marker={'opacity': 0.6, 'size': 10, 'color': 'orangered'},
                             showlegend=False,
                             ))
    traces.append(go.Scatter(x=non_sig_data['log2FC'],
                             y=non_sig_data['pval_plotting'],
                             text=non_sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Pass cutoff',
                             marker={'opacity': 0.6, 'size': 10, 'color': 'grey'},
                             showlegend=False
                             ))
    traces.append(go.Scatter(x=ticked_data['log2FC'],
                             y=ticked_data['pval_plotting'],
                             text=ticked_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers+text',
                             textposition='bottom center',
                             name='T',
                             marker={'opacity': 0.6, 'size': 10, 'color': 'green'},
                             showlegend=False
                             ))
    return {'data': traces,
            'layout': go.Layout(
                autosize=False,
                margin={
                    'l': 30,
                    'r': 15,
                    'pad': 0,
                    't': 30,
                    'b': 90, },
                # paper_bgcolor='rgba(0,0,0,0)',
                # plot_bgcolor = 'rgba(0,0,0,0)',
                xaxis={'title': 'log2FC'},
                yaxis={'title': '-log10(q-val)', 'ticktext': ticklab, 'tickvals': tickvals},
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
    return selected_data.to_dict('rows')


@app.callback(
    dash.dependencies.Output('dataset_table', 'data'),
    [dash.dependencies.Input('Sel_gene', 'value')])
def update_dataset_table(sel_gene):
    selected_data = main_data[main_data['Rv_ID']==sel_gene]
    selected_data['q-val'] = np.round(selected_data['q-val'], 2)
    selected_data['log2FC'] = np.round(selected_data['log2FC'], 2)
    selected_data = selected_data.sort_values(by='log2FC')
    return selected_data.to_dict('rows')


@app.callback(dash.dependencies.Output('cog', 'figure'),
              [dash.dependencies.Input('Sel_dataset', 'value'),
               dash.dependencies.Input('Sel_cog', 'value'),
               dash.dependencies.Input('log2FC', 'value'),
               dash.dependencies.Input('q-val', 'value')])
def update_cog(sel_dataset, sel_cog, log2FC, pval):
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    if sel_cog == 'Under-represented':
        sel_subset_filter = (selected_data['q-val'] <= pval) & (selected_data['log2FC'] <= -log2FC)
        colorscale = 'Cividis'
    else:
        sel_subset_filter = (selected_data['q-val'] <= pval) & (selected_data['log2FC'] >= -log2FC)
        colorscale = 'Viridis'
    sel_subset = selected_data[sel_subset_filter]
    cog_total_freq = cogs_df['COG'].value_counts(normalize=True)
    sel_cogs = cogs_df[cogs_df['X.Orf'].isin(sel_subset['Rv_ID'])]
    sel_cogs_freq = sel_cogs['COG'].value_counts(normalize=True)
    normalized_cogs = sel_cogs_freq / cog_total_freq
    normalized_cogs = normalized_cogs[~normalized_cogs.isnull()]
    normalized_cogs = normalized_cogs.sort_values()
    cog_names=cogs_desc.loc[list(normalized_cogs.index)]
    cog_names=list(cog_names.values)

    bar_data = [go.Bar(y=list(normalized_cogs.index), x=list(normalized_cogs.values),
                       orientation='h',
                       text=cog_names,
                       hoverinfo='text',
                       marker={'color': list(normalized_cogs.values), 'colorscale': colorscale})]
    return {'data': bar_data,
            'layout': go.Layout(
                margin={
                    'l': 30,
                    'r': 15,
                    'pad': 5,
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
    dash.dependencies.Output(component_id='orf_metadata', component_property='children'),
    [dash.dependencies.Input(component_id='Sel_gene', component_property='value')])
def print_orf_metadata(sel_gene):
    sel_details=orf_details[orf_details['ORF']==sel_gene]
    return sel_details['Description']

if __name__ == '__main__':
    app.run_server(debug=True)
