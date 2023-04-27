###随机生成一些点，然后可以动态显示出来他们的键。


import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
import numpy as np

# 创建一组随机点
num_points = 20
points = np.random.rand(num_points, 3)

# 创建Dash应用
app = dash.Dash(__name__)

# 定义应用布局
app.layout = html.Div([
    dcc.Graph(id='3d-scatter'),
    dcc.Slider(
        id='distance-slider',
        min=0,
        max=2,
        value=1,
        marks={i: f'{i}' for i in range(0, 3)},
        step=0.1
    )
])

# 计算点之间的距离
def distance(p1, p2):
    return np.sqrt(np.sum((p1 - p2) ** 2))

# 更新图形显示
@app.callback(
    Output('3d-scatter', 'figure'),
    [Input('distance-slider', 'value')])
def update_graph(max_distance):
    # 过滤出满足条件的连线
    lines = []
    for i in range(num_points):
        for j in range(i + 1, num_points):
            if distance(points[i], points[j]) <= max_distance:
                lines.append(go.Scatter3d(x=[points[i][0], points[j][0]],
                                          y=[points[i][1], points[j][1]],
                                          z=[points[i][2], points[j][2]],
                                          mode='lines',
                                          line=dict(color='rgba(0, 0, 255, 0.5)')))

    # 创建3D散点图
    scatter = go.Scatter3d(x=points[:, 0],
                           y=points[:, 1],
                           z=points[:, 2],
                           mode='markers',
                           marker=dict(size=6, color='red', opacity=0.8))

    # 更新图形数据和布局
    figure = go.Figure(data=[scatter] + lines)
    figure.update_layout(scene=dict(xaxis_title='X',
                                     yaxis_title='Y',
                                     zaxis_title='Z'))
    return figure

if __name__ == '__main__':
    app.run_server(debug=True)
