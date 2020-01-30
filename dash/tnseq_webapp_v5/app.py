import dash
import dash_core_components as dcc

external_scripts = ['https://code.jquery.com/jquery-3.3.1.min.js',
                    'https://cdn.datatables.net/v/dt/dt-1.10.18/datatables.min.js']
external_stylesheets = [dbc.themes.UNITED]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, external_scripts=external_scripts)
server=app.server
app.config.suppress_callback_exceptions = True
