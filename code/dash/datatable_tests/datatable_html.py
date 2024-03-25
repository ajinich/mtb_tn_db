import dash
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt


def Table(dataframe):
    rows = []
    for i in range(len(dataframe)):
        row = []
        for col in dataframe.columns:
            value = dataframe.iloc[i][col]
            # update this depending on which
            # columns you want to show links for
            # and what you want those links to be
            if col == 'Link':
                cell = html.Td(html.A(href='http://www.google.com', children='link_test'))
            else:
                cell = html.Td(children='link_test')
            row.append(cell)
        rows.append(html.Tr(row))
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        rows
    )






external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

test_table=pd.read_csv('for_column_test.csv')

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),
    dt.DataTable(
        id='table',
        columns=[{"name":x, "id":x} for x in test_table.columns],
        data=test_table.to_dict('records'),
        sort_action='native'
    ),
    Table(test_table),
    dcc.Graph(
        id='example-graph',
        figure={
            'data': [
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': 'SF'},
                {'x': [1, 2, 3], 'y': [2, 4, 5], 'type': 'bar', 'name': u'Montréal'},
            ],
            'layout': {
                'title': 'Dash Data Visualization'
            }
        }
    )

])

if __name__ == '__main__':
    app.run_server(debug=True)