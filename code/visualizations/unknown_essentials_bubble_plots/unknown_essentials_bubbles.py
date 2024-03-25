import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import os
from random import random
import seaborn as sns
import numpy as np
import plotly.graph_objs as go


def discretize_q_values(row, col_q_val):
    q_val = row[col_q_val]
    if q_val < 0.01:
        q_val_d = 3
    elif q_val < 0.05:
        q_val_d = 2
    else:
        q_val_d = 1
    return q_val_d


def unknown_essential_xy(TnSeq_screen, df_data, df_uk, rand_param = 0.6):

    # Grab data for a single TnSeq screen
    cols = [col for col in df_data.columns if TnSeq_screen in col]
    df_data_test = df_data[['Rv_ID', 'gene_name'] + cols].copy()
    
    # Discretize q-values: 
    col_q_val = [col for col in df_data_test.columns if 'q_val' in col][0]
    df_data_test['q_val_D'] = df_data_test.apply(discretize_q_values, 1, args=[col_q_val])
    
    # Merge with unknowns: 
    df_vis = df_data_test.merge(df_uk, on = ['Rv_ID', 'gene_name'], how = 'inner')
    
    # Get x-y datasets: 
    rv_ids = df_vis.Rv_ID.values
    uk_list = np.array(df_vis.UK_score_4)
    q_list = np.array(df_vis.q_val_D)
    
    # randomize: 
    uk_rd = np.array([uk + rand_param*random()-rand_param/2 for uk in uk_list])
    q_rd = np.array([q + rand_param*random()-rand_param/2 for q in q_list])
    
    # color the unknown-essentials differently: 
    current_palette = sns.color_palette()
    # all genes are gray by default. 
    color_list = np.array([(0.35,0.35,0.35)]*df_vis.shape[0])
    # Unknown essentials in a different color. 
    ind_temp = list(df_vis[(df_vis.q_val_D == 3) & (df_vis.UK_score_4 == 4)].index)
    color_list[ind_temp] = current_palette[0]
    color_list_rgb = ['rgb(' + ', '.join([str(np.round(rgb,2)) for rgb in col]) + ')' for col in color_list]
    
    return uk_rd, q_rd, color_list_rgb, rv_ids


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Load essentiality data:
path_data = '../../data/'
df_data = pd.read_excel(os.path.join(path_data, 'Tn_library_DB_qval_log2FC.xlsx'))
# Load annotation data: 
df_uk = pd.read_csv(os.path.join(path_data, 'unknown_essentials/unknown_ALL_levels_essential_scores.csv'))
df_uk = df_uk[['Rv_ID', 'gene_name', 'UK_score_4']]
# Grab essentiality / annotation data for a single TnSeq screen
TnSeq_screen = '2012_Zhang'
uk_rd, q_rd, color_list_rgb, rv_ids = unknown_essential_xy(TnSeq_screen, df_data, df_uk)

app.layout = html.Div([
    dcc.Graph(
        id='tSNE_plot',
        figure={
            'data': [
                go.Scatter(
                    x=uk_rd,
                    y=q_rd,
                    text=rv_ids,
                    mode='markers',
                    opacity=1,
                    hoverinfo = 'text',
                    marker={
                        'size': 15,
                        'line': {'width': 1.5, 'color': 'black'},
                        'color':color_list_rgb
                    }  
                ) 
            ],
            'layout': go.Layout(
                autosize=False,
    			width=800,
   				height=500,
    			xaxis = go.layout.XAxis(
                    tickmode = 'array',
                    tickvals = [0, 1, 2, 3, 4],
                    ticktext = ['most well\ncharacterized', '', '', '', 'least well\ncharacterized'],
                    tickfont=dict(size=14), 
                    title = 'Annotation'

                ),
                yaxis = go.layout.YAxis(
                    tickmode = 'array',
                    tickvals = [1, 2, 3],
                    ticktext = ['non-essential' ,'q-val < 0.05', 'q-val < 0.01'],
                    tickangle=270,
                    tickfont=dict(size=14),
                    title = 'Essentiality'
                ),
                margin={'l': 100, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)