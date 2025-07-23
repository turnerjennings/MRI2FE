import meshio

from ants import image_read

from ants.core.ants_image import ANTsImage

import os

import numpy as np

from ._MESHUTILS import mesh_wrapper


def nifti_to_inr(image: ANTsImage) -> str:
    # get dimensions and spacing
    xdim, ydim, zdim = image.shape

    vx, vy, vz = image.spacing

    # get data typeand length
    data = image.numpy()

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

    import tempfile

    temp_file = tempfile.NamedTemporaryFile(suffix=".inr", delete=False)

    out_path = temp_file.name

    try:
        temp_file.write(header.encode("ascii"))
        temp_file.write(data.tobytes(order="F"))
        temp_file.close()

        if not os.path.exists(out_path):
            raise ValueError(f"File was not written to {out_path}")

        return out_path

    except:
        if os.path.exists(out_path):
            os.unlink(out_path)
        raise


def mesh_from_nifti(
    filepath: str,
    optimize: bool = False,
    facetAngle: float = 30.0,
    facetSize: float = 1.0,
    facetDistance: float = 4.0,
    cellRadiusEdgeRatio: float = 3.0,
    cellSize: float = 1.0,
) -> meshio.Mesh:
    if not os.path.exists(filepath):
        raise ValueError(f"input filepath {filepath} does not exist")

    image = image_read(filepath)

    origin = image.origin

    transfer_path = nifti_to_inr(image=image)

    mesh_path = mesh_wrapper(
        transfer_path,
        optimize,
        facetAngle,
        facetSize,
        facetDistance,
        cellRadiusEdgeRatio,
        cellSize,
    )

    if not os.path.exists(mesh_path):
        raise ValueError(f"Mesh wrapper did not create mesh at {mesh_path}")

    try:
        mesh: meshio.Mesh = meshio.read(mesh_path)

        mesh.points = mesh.points - np.array(origin)

        return mesh
    except ValueError:
        raise ValueError(
            "Error loading mesh from file, mesh may be too large..."
        )
