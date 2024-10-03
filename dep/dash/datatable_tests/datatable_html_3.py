import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import visdcc
import pandas as pd

test_table=pd.read_csv('for_column_test.csv')

def Table(dataframe, id):
    rows = []
    for i in range(len(dataframe)):
        row = []
        for col in dataframe.columns:
            value = dataframe.iloc[i][col]
            # update this depending on which
            # columns you want to show links for
            # and what you want those links to be
            if col == 'Link':
                cell = html.Td(html.A(href='http://www.google.com', children=value))
            else:
                cell = html.Td(children=value)
            row.append(cell)
        rows.append(html.Tr(row))
    Thead = html.Thead(
        [html.Tr([html.Th(col, style={'background-color':"#FF0000"}) for col in dataframe.columns])]
    )
    Tbody=html.Tbody(rows)

    return html.Table([Thead, Tbody], id=id, className="display")



external_scripts = ['https://code.jquery.com/jquery-3.3.1.min.js',
                    'https://cdn.datatables.net/v/dt/dt-1.10.18/datatables.min.js']
external_stylesheets = ['https://cdn.datatables.net/v/dt/dt-1.10.18/datatables.min.css']

app = dash.Dash(external_scripts=external_scripts,
                external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Dropdown(
        id='demo-dropdown',
        options=[
            {'label': 'Anisha', 'value': 'Anisha'},
            {'label': 'Adrian', 'value': 'Adrian'},
        ],
        value='Anisha'
    ),
    visdcc.Run_js(id='javascript', run="$('#datatable1').DataTable()"),
    html.Div(id='output_div')
])

@app.callback(
    Output('output_div', 'children'),
    [Input('demo-dropdown', 'value')])
def return_table(selected_value):
    df=test_table[test_table['Name']==selected_value]
    return Table(df, id='datatable1')


# @app.callback(
#     Output('javascript', 'run'),
#     [Input('button', 'n_clicks')])
# def myfun(x):
#     return "$('#datatable').DataTable()"


# @app.callback(
#     Output('output_div', 'children'),
#     [Input('datatable_{}_{}'.format(i, j), 'n_clicks') for i in range(len(df)) for j in range(len(df.columns))])
# def myfun(*args):
#     ctx = dash.callback_context
#     if not ctx.triggered or ctx.triggered[0]['value'] is None:  return ""
#     input_id = ctx.triggered[0]['prop_id'].split('.')[0]
#     row = input_id.split('_')[1]
#     col = input_id.split('_')[2]
#     return "you click on row {} col {}".format(row, col)
#

if __name__ == '__main__':
    app.run_server(debug=True)