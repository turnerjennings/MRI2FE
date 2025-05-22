import meshio

from _MESHUTILS import mesh_wrapper



def mesh_from_nifti(
        filepath:str,
        optimize:bool = False,
        facetAngle:float = 30.,
        facetSize:float = 1.,
        facetDistance:float = 4.,
        cellRadiusEdgeRatio:float = 3.,
        cellSize:float = 1.
) -> meshio.Mesh:
    
    out_path = mesh_wrapper(
        filepath,
        optimize,
        facetAngle,
        facetSize,
        facetDistance,
        cellRadiusEdgeRatio,
        cellSize
    )

    mesh = meshio.read(out_path)

    return mesh