import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def icosahedron_vertices():
    phi = (1 + np.sqrt(5)) / 2
    vertices = np.array([
        [-1, phi, 0],
        [1, phi, 0],
        [-1, -phi, 0],
        [1, -phi, 0],
        [0, -1, phi],
        [0, 1, phi],
        [0, -1, -phi],
        [0, 1, -phi],
        [phi, 0, -1],
        [phi, 0, 1],
        [-phi, 0, -1],
        [-phi, 0, 1]
    ], dtype=float)
    return vertices

def icosahedron_faces():
    # Faces (vertex indices)
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)
    ]

    return faces

def subdivide_triangle(p1, p2, p3, depth):
    if depth == 0:
        return [p1, p2, p3]
    
    # Midpoints of the edges
    m1 = (p1 + p2) / 2
    m2 = (p2 + p3) / 2
    m3 = (p3 + p1) / 2

    # New vertices
    v1 = p1
    v2 = p2
    v3 = p3

    # Recursively subdivide
    triangles = []
    triangles.extend(subdivide_triangle(v1, m1, m3, depth - 1))
    triangles.extend(subdivide_triangle(m1, v2, m2, depth - 1))
    triangles.extend(subdivide_triangle(m3, m2, v3, depth - 1))
    triangles.extend(subdivide_triangle(m1, m2, m3, depth - 1))

    return triangles

def generate_icosahedral_grid(depth):
    icosahedron_vertices_array = icosahedron_vertices()
    faces = icosahedron_faces()
    triangles = []
    
    # Subdivide each face of the icosahedron
    for face in faces:
        p1, p2, p3 = icosahedron_vertices_array[face[0]], icosahedron_vertices_array[face[1]], icosahedron_vertices_array[face[2]]
        triangles.extend(subdivide_triangle(p1, p2, p3, depth))

    return np.array(triangles)

# def project_onto_surface(vertices):
    
    

def plot_icosahedral_grid(vertices):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Icosahedral Grid')

    plt.show()

if __name__ == "__main__":
    depth = 4  # Adjust the depth of subdivision
    icosahedral_grid = generate_icosahedral_grid(depth)
    icosahedral_grid = np.unique(icosahedral_grid, axis=0)
    # project_onto_surface(icosahedral_grid)
    for i in range(len(icosahedral_grid)):
        icosahedral_grid[i] = normalize(icosahedral_grid[i])
    
    for i, v in enumerate(icosahedral_grid):
        if np.linalg.norm(v - np.array([0,0,1])) <1e-8:
            icosahedral_grid = np.delete(icosahedral_grid, i, axis=0)
            break
    for i, v in enumerate(icosahedral_grid):
        if np.linalg.norm(v - np.array([0,0,-1])) <1e-8:
            icosahedral_grid = np.delete(icosahedral_grid, i, axis=0)
            break

    print(len(icosahedral_grid))
    np.savetxt(f'icosa_d{depth}_N{len(icosahedral_grid)}.dat', icosahedral_grid)
    # plot_icosahedral_grid(icosahedral_grid)
