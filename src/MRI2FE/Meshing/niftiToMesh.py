from ants.core.ants_image import ANTsImage

from skimage.measure import marching_cubes

import numpy as np

import tetgen

import pyvista as pv

def find_surface(img:ANTsImage, 
                 low_value:float, 
                 high_value:float)->pv.PolyData:

    img_data = img.numpy()

    img_spacing = img.spacing

    if len(img_data.shape) > 3:
        raise ValueError("Image must be 3D")
    
    threshold = 0.5*(high_value - low_value) + low_value
    
    verts, faces, _, _ = marching_cubes(volume=img_data,
                                        level=threshold,
                                        spacing=img_spacing)
    
    faces_pv = np.insert(faces, 0, 3, axis=1).flatten()
    
    surf = pv.PolyData(verts,faces_pv)

    if surf.is_manifold is False:
        raise ValueError("Generated surface is not manifold")

    return surf

def mesh_between_surface(mesh_o:pv.PolyData,
                         mesh_i:pv.PolyData = None) -> np.ndarray:
    
    
    
    #combine inside and outside surfaces

    if not mesh_o.is_manifold:
        raise ValueError("Outside surface is not manifold")
    
    if mesh_i is not None and not mesh_i.is_manifold:
        raise ValueError("Inside surface is not manifold")
    
    if mesh_i is None:
        tgen = tetgen.TetGen(mesh_o)
        tgen.make_manifold()
        nodes, elements = tgen.tetrahedralize()
        return nodes, elements
    
    mesh_i.flip_faces(inplace=True)
    
    def extract_triangles(mesh):
        """Extract triangular faces from PyVista mesh."""
        vtk_faces = mesh.faces
        triangles = []
        i = 0
        while i < len(vtk_faces):
            n_points = vtk_faces[i]
            if n_points == 3:  # We're only dealing with triangles
                triangles.append([vtk_faces[i+1], vtk_faces[i+2], vtk_faces[i+3]])
            i += n_points + 1
        return np.array(triangles)
    
    # Extract triangular faces
    faces_o = extract_triangles(mesh_o)
    faces_i = extract_triangles(mesh_i)
    
    # If the inner mesh is meant to be a hole, flip its normals
    # This is important for defining the inside/outside of the domain
    faces_i = faces_i[:, [0, 2, 1]]  # Flip winding order
    
    # Combine points and adjust inner face indices
    n_outer_points = mesh_o.n_points
    adjusted_faces_i = faces_i.copy()
    adjusted_faces_i += n_outer_points  # Shift indices for inner faces
    
    combined_points = np.vstack([mesh_o.points, mesh_i.points])
    
    # Format for PyVista - combine the faces
    # First, reformat triangles for PyVista
    outer_faces_pv = np.hstack([np.ones((faces_o.shape[0], 1), dtype=int) * 3, faces_o])
    inner_faces_pv = np.hstack([np.ones((adjusted_faces_i.shape[0], 1), dtype=int) * 3, adjusted_faces_i])
    
    # Combine and flatten
    combined_faces = np.vstack([outer_faces_pv, inner_faces_pv]).flatten()
    
    # Create the combined mesh
    combined_surf = pv.PolyData(combined_points, combined_faces)
    
    # Check if combined surface is still manifold
    if not combined_surf.is_manifold:
        print("Warning: Combined surface is not manifold. Attempting to fix...")
        try:
            import pymeshfix
            meshfix = pymeshfix.MeshFix(combined_surf)
            meshfix.repair()
            combined_surf = meshfix.mesh
            if not combined_surf.is_manifold:
                raise ValueError("Failed to make combined surface manifold")
        except ImportError:
            print("pymeshfix not available to repair mesh. Trying TetGen's built-in repair...")
    
    # Create TetGen instance and generate the tetrahedral mesh
    tgen = tetgen.TetGen(combined_surf)

    try:
        # Try to generate the tetrahedral mesh
        nodes, elements = tgen.tetrahedralize(
            order=1,          # Linear tetrahedra
            mindihedral=10,   # Minimum dihedral angle
            minratio=1.5,     # Quality constraint
            quality=True      # Enable quality mesh generation
        )
        return nodes, elements
    except RuntimeError as e:
        # Provide more diagnostic information
        print(f"TetGen failed: {e}")
        print("Mesh statistics:")
        print(f"  Outer mesh: {mesh_o.n_points} points, {mesh_o.n_faces} faces")
        if mesh_i is not None:
            print(f"  Inner mesh: {mesh_i.n_points} points, {mesh_i.n_faces} faces")
        print(f"  Combined mesh: {combined_surf.n_points} points, {combined_surf.n_faces} faces")
        
        # Visualize to help diagnose
        p = pv.Plotter()
        p.add_mesh(mesh_o, color='blue', opacity=0.5, label='Outer')
        if mesh_i is not None:
            p.add_mesh(mesh_i, color='red', opacity=0.5, label='Inner')
        p.add_legend()
        p.show()
        
        # Re-raise with additional context
        raise ValueError(f"Failed to tetrahedralize: {e}")
