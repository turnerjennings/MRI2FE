import pytest

from ants import from_numpy

import numpy as np

from MRI2FE.Meshing import find_surface,mesh_between_surface

import os

import pyvista as pv

def plot_tetrahedral_mesh(fpath:str, nodes, elements, outer_surface=None, inner_surface=None, 
                          show_edges=True):
    """
    Plot a tetrahedral mesh and optionally its boundary surfaces.
    
    Parameters:
    -----------
    nodes : np.ndarray
        Array of shape (n_nodes, 3) containing the node coordinates
    elements : np.ndarray
        Array of shape (n_elements, 4) containing the tetrahedron connectivity
    outer_surface : pv.PolyData, optional
        The outer boundary surface
    inner_surface : pv.PolyData, optional
        The inner boundary surface
    show_edges : bool, default=True
        Whether to show edges of the tetrahedral elements
    opacity : float, default=0.7
        Opacity of the tetrahedral mesh (0-1)
    cmap : str, default='viridis'
        Colormap for the tetrahedral mesh
    """
    import pyvista as pv
    
    # Create the tetrahedral mesh as an UnstructuredGrid
    # Convert elements to the format expected by PyVista if needed
    # PyVista expects VTK_TETRA (10) as the cell type for tetrahedra
    if elements.shape[1] == 4:  # If elements are just the 4 node indices
        vtk_elements = np.hstack([
            np.ones((elements.shape[0], 1), dtype=np.int64) * 4,  # 4 points per tetrahedron
            elements
        ])
        # Flatten the array as PyVista expects a 1D array
        vtk_elements = vtk_elements.flatten()
    else:
        vtk_elements = elements  # Assume already in the correct format
    
    # Create the unstructured grid
    tet_mesh = pv.UnstructuredGrid(vtk_elements, [10] * (elements.shape[0]), nodes)
    
    # Compute a quality metric for the cells (optional)
    tet_mesh.compute_cell_quality('collapse_ratio')
    
    tet_mesh_clipped = tet_mesh.clip('z')

    # Create a plotter
    p = pv.Plotter()
    
  
    p.add_mesh(tet_mesh_clipped,
               show_edges=show_edges, 
               scalar_bar_args={'title': 'Cell Quality'})
    
    # Add boundary surfaces if provided
    if outer_surface is not None:
        p.add_mesh(outer_surface, color='blue', opacity=0.3, label='Outer Surface')
    
    if inner_surface is not None:
        p.add_mesh(inner_surface, color='red', opacity=0.3, label='Inner Surface')
    
    # Add a legend if surfaces are shown
    if outer_surface is not None or inner_surface is not None:
        p.add_legend()

    p.save_graphic(fpath)

def generate_MRI_sphere(radius:float):

    array_shape = (101,101,101)

    spacing = (1.,1.,1.)

    origin = (0.,0.,0.)

    direction = np.eye(3)

    arr = np.zeros(array_shape)


    center = np.array([50, 50, 50])

    # Create coordinate grids
    z, y, x = np.ogrid[:array_shape[0], :array_shape[1], :array_shape[2]]
    
    # Calculate squared distance from center for each point
    # (Using squared distance avoids computing square roots)
    dist_squared = (z - center[0])**2 + (y - center[1])**2 + (x - center[2])**2
    
    # Set points within the radius to 1
    mask = dist_squared <= radius**2
    arr[mask] = 1

    sphere=from_numpy(arr,origin=origin, spacing=spacing,direction=direction)

    return sphere

class TestNiftiToMesh:

    def test_sphere_surface(self):
        image = generate_MRI_sphere(20)

        mesh = find_surface(img=image, low_value=0., high_value=1.)

        print(mesh)

        for i in range(mesh.verts.shape[0]):
            dist = np.linalg.norm(mesh.verts[i,:]-np.array([50.,50.,50.]))
            assert (dist - 20)**2 <= 0.5
    

    def test_mesh_solid(self):
        image = generate_MRI_sphere(20.0)

        mesh = find_surface(img = image,
                                    low_value=0.,
                                    high_value=1.)
        nodes,els = mesh_between_surface(mesh_o = mesh)

        for i in range(nodes.shape[0]):
            dist = np.linalg.norm(nodes[i,:]-np.array([50.,50.,50.]))
            assert dist <= 20.5


    def test_mesh_between(self):
        image1 = generate_MRI_sphere(20.0)
        image2 = generate_MRI_sphere(5.0)

        mesh_o = find_surface(img = image1,
                                    low_value=0.,
                                    high_value=1.)
        
        mesh_i = find_surface(img = image2,
                                    low_value=0.,
                                    high_value=1.)
        
        for i in range(mesh_o.points.shape[0]):
            dist = np.linalg.norm(mesh_o.points[i,:]
                                   - np.array([50.,50.,50.]))
            assert (dist - 20.)**2 <= 0.5

        for i in range(mesh_i.points.shape[0]):
            dist = np.linalg.norm(mesh_i.points[i,:]
                                   - np.array([50.,50.,50.]))
            assert (dist - 5.)**2 <= 0.5
        
        nodes, els = mesh_between_surface(mesh_o,mesh_i)
        #create pyvista mesh



        #plot resulting mesh
        if not os.path.exists("./.benchmark/"):
            os.mkdir("./.benchmark/")

        plot_tetrahedral_mesh("./.benchmark/test_mesh_between.pdf",
                              nodes,
                              els)

        for i in range(nodes.shape[0]):
            dist = np.linalg.norm(nodes[i,:]
                                   - np.array([50.,50.,50.]))
            
            assert dist <= 20.5
            assert dist >= 4.5

        
        
        






