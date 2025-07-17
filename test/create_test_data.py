import numpy as np
from ants import from_numpy, image_write
import os
from ..src.MRI2FE.generate_mesh import mesh_from_nifti
import meshio


def create_test_nifti_MRI_files(output_dir: str = "test/test_data"):
    """Create tets NIfTI file for geometry and the associated .inr file

    Args:
        output_dir (str, optional): _description_. Defaults to "test/test_data".
    """

    # Parameters
    shape = (64, 64, 64)
    center = np.array(shape) // 2
    radii = [5, 10, 15, 20]  # Radii for spheres
    labels = list(range(1, len(radii) + 1))  # Labels: 1, 2, 3, ...

    spacing = (1.0, 1.0, 1.0)
    origin = (4.0, -1.0, 2.0)

    # Create empty volume
    volume = np.zeros(shape, dtype=np.uint8)

    # Grid of coordinates
    zz, yy, xx = np.meshgrid(
        np.arange(shape[0]),
        np.arange(shape[1]),
        np.arange(shape[2]),
        indexing="ij",
    )
    dist = np.sqrt(
        (xx - center[2]) ** 2 + (yy - center[1]) ** 2 + (zz - center[0]) ** 2
    )

    # Assign labels to voxels inside each sphere
    for i, r in enumerate(radii):
        if i == 0:
            mask = dist <= r
        else:
            mask = (dist <= r) & (dist > radii[i - 1])
        volume[mask] = labels[i]

    # Save as NIfTI
    nifti_img = from_numpy(data=volume, origin=origin, spacing=spacing)

    image_write(
        nifti_img, os.path.join(output_dir, "test_concentric_spheres.nii")
    )

    # get dimensions and spacing
    xdim, ydim, zdim = nifti_img.shape

    vx, vy, vz = nifti_img.spacing

    # get data typeand length

    data = nifti_img.numpy()

    btype, bitlen = {
        "uint8": ("unsigned fixed", 8),
        "uint16": ("unsigned fixed", 16),
        "float32": ("float", 32),
        "float64": ("float", 64),
    }[data.dtype.name]

    # encode file
    header = "\n".join(
        [
            "#INRIMAGE-4#{",
            f"XDIM={xdim}",
            f"YDIM={ydim}",
            f"ZDIM={zdim}",
            "VDIM=1",
            f"TYPE={btype}",
            f"PIXSIZE={bitlen} bits",
            "CPU=decm",
            f"VX={vx:f}",
            f"VY={vy:f}",
            f"VZ={vz:f}",
        ]
    )
    header += "\n"

    header = header + "\n" * (256 - 4 - len(header)) + "##}\n"

    out_path = os.path.join(output_dir, "test_concentric_spheres.inr")

    with open(out_path, "wb") as temp_file:  # Note: 'wb' instead of 'w'
        temp_file.write(header.encode("ascii"))
        temp_file.write(data.tobytes(order="F"))

    # generate mesh
    mesh_obj: meshio.Mesh = mesh_from_nifti(
        os.path.join(output_dir, "test_concentric_spheres.nii")
    )
    mesh_obj.write(os.path.join(output_dir, "test_mesh_out.mesh"))


def create_test_nifti_MRE_files(output_dir: str = "test/test_data"):
    """Create test NIfTI files for stiffness and damping ratio.

    Args:
        output_dir: Directory to save the test files
    """
    os.makedirs(output_dir, exist_ok=True)
    shape = (20, 20, 20)
    stiffness = np.ones(shape) * 5.0
    x, y, z = np.meshgrid(
        np.linspace(-1, 1, shape[0]),
        np.linspace(-1, 1, shape[1]),
        np.linspace(-1, 1, shape[2]),
    )
    r = np.sqrt(x**2 + y**2 + z**2)
    stiffness += 3.0 * np.exp(-(r**2))

    damping = np.ones(shape) * 0.2
    damping += 0.1 * np.exp(-(r**2))

    stiffness_img = from_numpy(stiffness)
    damping_img = from_numpy(damping)

    image_write(stiffness_img, os.path.join(output_dir, "test_stiffness.nii"))
    image_write(
        damping_img, os.path.join(output_dir, "test_damping_ratio.nii")
    )

    return os.path.join(output_dir, "test_stiffness.nii"), os.path.join(
        output_dir, "test_damping_ratio.nii"
    )


def create_test_mesh_file(output_dir: str = "test/test_data"):
    """Create a simple tetrahedral mesh file for testing.

    Args:
        output_dir: Directory to save the test files
    """
    os.makedirs(output_dir, exist_ok=True)

    nodes = np.array(
        [
            [1, 0, 0, 0],
            [2, 1, 0, 0],
            [3, 0, 1, 0],
            [4, 1, 1, 0],
            [5, 0, 0, 1],
            [6, 1, 0, 1],
            [7, 0, 1, 1],
            [8, 1, 1, 1],
        ]
    )

    elements = np.array(
        [
            [1, 4, 1, 2, 3, 5, 0, 0, 0, 0, 0, 0],
            [2, 4, 2, 4, 3, 6, 0, 0, 0, 0, 0, 0],
            [3, 4, 3, 4, 7, 5, 0, 0, 0, 0, 0, 0],
            [4, 5, 4, 8, 7, 6, 0, 0, 0, 0, 0, 0],
            [5, 5, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
        ]
    )

    output_file = os.path.join(output_dir, "test_mesh.k")
    with open(output_file, "w") as f:
        f.write("*KEYWORD\n")
        f.write("*TITLE\nTest Mesh\n\n")

        f.write("*NODE\n")
        f.write(
            "$#   nid               x               y               z      tc      rc\n"
        )
        for node in nodes:
            f.write(
                f"{int(node[0]):8d}{node[1]:16.6f}{node[2]:16.6f}{node[3]:16.6f}\n"
            )

        f.write("\n*ELEMENT_SOLID\n")
        f.write(
            "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n"
        )
        for elem in elements:
            f.write(
                f"{int(elem[0]):8d}{int(elem[1]):8d}{int(elem[2]):8d}{int(elem[3]):8d}{int(elem[4]):8d}{int(elem[5]):8d}{int(elem[6]):8d}{int(elem[7]):8d}{int(elem[8]):8d}{int(elem[9]):8d}\n"
            )

        f.write("\n*PART\n")
        f.write(
            "$#                                                                         title\n"
        )
        f.write(
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n"
        )
        f.write(f"{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{0:9d}\n")

        f.write("\n*PART\n")
        f.write(
            "$#                                                                         title\n"
        )
        f.write(
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n"
        )
        f.write(f"{2:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{0:9d}\n")

        f.write("\n*SECTION_SOLID\n")
        f.write("$#   secid    elform       aet\n")
        f.write("         1         1         0\n")

        f.write("\n*END\n")

    return output_file


def create_meshio_mesh(output_dir: str = "test/test_data"):
    path = os.path.join(output_dir, "test_concentric_spheres.nii")

    mesh = mesh_from_nifti(path, optimize=True)

    out_path = os.path.join(output_dir, "test_mesh_out.mesh")

    mesh.write(out_path)

    return out_path


def create_all_test_files(output_dir: str = "test/test_data"):
    """Create all test files needed for the MRE pipeline.

    Args:
        output_dir: Directory to save the test files
    """
    stiffness_path, damping_path = create_test_nifti_MRE_files(output_dir)
    mesh_path = create_test_mesh_file(output_dir)
    meshio_path = create_meshio_mesh(output_dir)

    print(f"Created test files:")
    print(f"Stiffness file: {stiffness_path}")
    print(f"Damping ratio file: {damping_path}")
    print(f"Keyword Mesh file: {mesh_path}")
    print(f"Mesh file: {meshio_path}")

    return stiffness_path, damping_path, mesh_path


if __name__ == "__main__":
    create_all_test_files()
