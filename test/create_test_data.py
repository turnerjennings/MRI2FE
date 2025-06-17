import numpy as np
from ants import from_numpy, image_write
import os

def create_test_nifti_files(output_dir: str = "test/test_data"):
    """Create test NIfTI files for stiffness and damping ratio.
    
    Args:
        output_dir: Directory to save the test files
    """
    os.makedirs(output_dir, exist_ok=True)
    shape = (20, 20, 20)
    stiffness = np.ones(shape) * 5.0
    x, y, z = np.meshgrid(np.linspace(-1, 1, shape[0]),
                         np.linspace(-1, 1, shape[1]),
                         np.linspace(-1, 1, shape[2]))
    r = np.sqrt(x**2 + y**2 + z**2)
    stiffness += 3.0 * np.exp(-r**2)
    
    damping = np.ones(shape) * 0.2
    damping += 0.1 * np.exp(-r**2)
    
    stiffness_img = from_numpy(stiffness)
    damping_img = from_numpy(damping)
    
    image_write(stiffness_img, os.path.join(output_dir, "test_stiffness.nii"))
    image_write(damping_img, os.path.join(output_dir, "test_damping_ratio.nii"))
    
    return os.path.join(output_dir, "test_stiffness.nii"), os.path.join(output_dir, "test_damping_ratio.nii")

def create_test_mesh_file(output_dir: str = "test/test_data"):
    """Create a simple tetrahedral mesh file for testing.
    
    Args:
        output_dir: Directory to save the test files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    nodes = np.array([
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 0, 1, 0],
        [4, 1, 1, 0],
        [5, 0, 0, 1],
        [6, 1, 0, 1],
        [7, 0, 1, 1],
        [8, 1, 1, 1],
    ])
    
    elements = np.array([
        [1, 4, 1, 2, 3, 5, 0, 0, 0, 0, 0, 0],
        [2, 4, 2, 4, 3, 6, 0, 0, 0, 0, 0, 0],
        [3, 4, 3, 4, 7, 5, 0, 0, 0, 0, 0, 0],
        [4, 5, 4, 8, 7, 6, 0, 0, 0, 0, 0, 0],
        [5, 5, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    ])
    
    output_file = os.path.join(output_dir, "test_mesh.k")
    with open(output_file, 'w') as f:
        f.write("*KEYWORD\n")
        f.write("*TITLE\nTest Mesh\n\n")
        
        f.write("*NODE\n")
        f.write("$#   nid               x               y               z      tc      rc\n")
        for node in nodes:
            f.write(f"{int(node[0]):8d}{node[1]:16.6f}{node[2]:16.6f}{node[3]:16.6f}\n")
        
        f.write("\n*ELEMENT_SOLID\n")
        f.write("$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n")
        for elem in elements:
            f.write(f"{int(elem[0]):8d}{int(elem[1]):8d}{int(elem[2]):8d}{int(elem[3]):8d}{int(elem[4]):8d}{int(elem[5]):8d}{int(elem[6]):8d}{int(elem[7]):8d}{int(elem[8]):8d}{int(elem[9]):8d}\n")
        
        f.write("\n*PART\n")
        f.write("$#                                                                         title\n")
        f.write("$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n")
        f.write(f"{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{0:9d}\n")
        
        f.write("\n*PART\n")
        f.write("$#                                                                         title\n")
        f.write("$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n")
        f.write(f"{2:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{0:9d}\n")
        
        f.write("\n*SECTION_SOLID\n")
        f.write("$#   secid    elform       aet\n")
        f.write("         1         1         0\n")
        
        f.write("\n*END\n")
    
    return output_file

def create_all_test_files(output_dir: str = "test/test_data"):
    """Create all test files needed for the MRE pipeline.
    
    Args:
        output_dir: Directory to save the test files
    """
    stiffness_path, damping_path = create_test_nifti_files(output_dir)
    mesh_path = create_test_mesh_file(output_dir)
    
    print(f"Created test files:")
    print(f"Stiffness file: {stiffness_path}")
    print(f"Damping ratio file: {damping_path}")
    print(f"Mesh file: {mesh_path}")
    
    return stiffness_path, damping_path, mesh_path

if __name__ == "__main__":
    create_all_test_files() 