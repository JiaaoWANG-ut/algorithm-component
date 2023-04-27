import plotly.graph_objects as go
import numpy as np
import plotly.offline as offline

def plot_lattice(A, rand_coords):
    # Define the vertices and edges of the unit cell
    x = [0, 1, 1, 0, 0, 1, 1, 0]
    y = [0, 0, 1, 1, 0, 0, 1, 1]
    z = [0, 0, 0, 0, 1, 1, 1, 1]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]

    # Generate 10 random points inside the unit cell
    num_points = 10
    #rand_coords = np.random.rand(num_points, 3)

    # Create the figure
    fig = go.Figure()

    # Add the unit cell and its edges
    fig.add_trace(go.Mesh3d(x=x, y=y, z=z, color='lightblue', opacity=0.5))
    for edge in edges:
        x0, y0, z0 = A @ np.array([x[edge[0]], y[edge[0]], z[edge[0]]])
        x1, y1, z1 = A @ np.array([x[edge[1]], y[edge[1]], z[edge[1]]])
        fig.add_trace(go.Scatter3d(x=[x0, x1], y=[y0, y1], z=[z0, z1], mode='lines', line=dict(color='blue', width=5)))

    # Add the random points inside the unit cell
    for point in rand_coords:
        transformed_point = A @ point
        fig.add_trace(go.Scatter3d(x=[transformed_point[0]], y=[transformed_point[1]], z=[transformed_point[2]], 
                                   mode='markers', marker=dict(size=6, color='red')))

    # Save the visualization as an HTML file
    offline.plot(fig, filename='lattice.html', auto_open=True)


#A = np.array([[1, 3, 0], [0, 1, 3], [0, 0, 1]])
#rand_coords = np.random.rand(10, 3)*2
#plot_lattice(A, rand_coords)
