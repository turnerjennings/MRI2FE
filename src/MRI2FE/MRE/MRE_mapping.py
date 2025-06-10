import numpy as np
from ants.core.ants_image import ANTsImage
import scipy.spatial as sp
from ..utilities import COM_align, spatial_map
from ..FEModel.femodel import FEModel
import os
from datetime import datetime


def map_MRE_to_mesh(
    fe_model: FEModel,
    map: ANTsImage,
    elcentroids: np.ndarray,
    ect: np.ndarray,
    offset: int = 3,
    csvpath: str = None,
):
    """Map elements to parts from a segmented spatial map

    Args:
        fe_model (FEModel): Finite element model object
        map (AntsImage): segmented spatial map in voxel space
        elcentroids (np.ndarray): (3,n_elements) array of the coordinates of each element centroid
        ect (np.ndarray): element connectivity table in 10-node LS-Dyna format
        offset (int, optional): optional filter if it is not desired to remap all pids, will skip mapping any elements belonging to pid <= offset. Defaults to 3.
        csvpath (str, optional): Path to save CSV output files. Defaults to None.

    Raises:
        TypeError: If input types are invalid
        ValueError: If input dimensions or values are invalid
        FileNotFoundError: If CSV directory doesn't exist

    Returns:
        ect (np.ndarray): new ect in 10-node LS-DYNA format with updated PIDs
    """
    if not isinstance(fe_model, FEModel):
        raise TypeError("fe_model must be a FEModel object")

    if not isinstance(map, ANTsImage):
        raise TypeError("map must be an ANTsImage object")

    if not isinstance(elcentroids, np.ndarray):
        raise TypeError("elcentroids must be a numpy array")
    if not isinstance(ect, np.ndarray):
        raise TypeError("ect must be a numpy array")

    if len(elcentroids.shape) != 2 or elcentroids.shape[1] != 3:
        raise ValueError(
            "elcentroids must be a 2D array with shape (n_elements, 3)"
        )
    if len(ect.shape) != 2:
        raise ValueError("ect must be a 2D array")

    if not isinstance(offset, int):
        raise TypeError("offset must be an integer")
    if offset < 0:
        raise ValueError("offset must be non-negative")

    if csvpath is not None:
        if not isinstance(csvpath, str):
            raise TypeError("csvpath must be a string")
        csv_dir = os.path.dirname(csvpath)
        if csv_dir and not os.path.exists(csv_dir):
            raise FileNotFoundError(
                f"CSV output directory does not exist: {csv_dir}"
            )

    physical_space_map = spatial_map(map)
    # print(physical_space_map)

    # create filtered space map and centroids for brain only
    physical_space_map_nonzero = physical_space_map[
        physical_space_map[:, 3] > 0
    ]
    physical_space_map_nonzero[:, [0, 1]] = physical_space_map_nonzero[
        :, [1, 0]
    ]
    # print(
    #    f"physical space map: min/max x=({np.min(physical_space_map_nonzero[:,0])},{np.max(physical_space_map_nonzero[:,0])}), y=({np.min(physical_space_map_nonzero[:,1])},{np.max(physical_space_map_nonzero[:,1])}), z=({np.min(physical_space_map_nonzero[:,2])},{np.max(physical_space_map_nonzero[:,2])})"
    # )

    ect_brain_mask = ect[:, 1] > 3
    elcentroids_brain = elcentroids[ect_brain_mask, :]

    # print(
    #    f"physical_space_map_nonzero shape: {physical_space_map_nonzero.shape}\nelcentroids_brain shape: {elcentroids_brain.shape}"
    # )
    print(f"{datetime.now()}\t\tAligning COM...")
    physical_space_map_transformed = COM_align(
        elcentroids,
        physical_space_map_nonzero[:, 0:3],
        fixed_mask=elcentroids_brain,
    )

    physical_space_map_nonzero[:, 0:3] = physical_space_map_transformed

    # print(
    #    f"Transformed physical space map: min/max x=({np.min(physical_space_map_nonzero[:,0])},{np.max(physical_space_map_nonzero[:,0])}), y=({np.min(physical_space_map_nonzero[:,1])},{np.max(physical_space_map_nonzero[:,1])}), z=({np.min(physical_space_map_nonzero[:,2])},{np.max(physical_space_map_nonzero[:,2])})"
    # )
    # write csv output if requested
    if csvpath is not None:
        elcentroids_brain_out = np.hstack(
            (
                elcentroids_brain,
                (np.max(physical_space_map_nonzero[:, 3]) + 1)
                * np.ones((elcentroids_brain.shape[0], 1)),
            )
        )
        np.savetxt(
            csvpath + "_braincentroids.csv",
            elcentroids_brain_out,
            delimiter=",",
            header="X,Y,Z,ID",
        )
        np.savetxt(
            csvpath + "_MREvoxels.csv",
            physical_space_map_nonzero,
            delimiter=",",
            header="X,Y,Z,ID",
        )

    # create KDTree
    physical_space_tree = sp.KDTree(
        physical_space_map_nonzero[:, 0:3], leafsize=15
    )

    # print(f"physical space map: number of points = {physical_space_tree.n}, number of dimensions = {physical_space_tree.m}, number of nodes = {physical_space_tree.size}")

    query_mask = ect[:, 1] > offset

    element_query = elcentroids[query_mask, :]

    d, idx = physical_space_tree.query(element_query, k=1)

    # idx_avg = np.round(np.mean(idx,axis=1))

    new_pids = physical_space_map_nonzero[idx, 3] + offset
    # np.savetxt("test_query.csv",np.hstack((new_pids,d)),delimiter=',')

    ect[query_mask, 1] = new_pids

    # Update element connectivity table in FEModel
    for element_id, part_id in enumerate(ect[:, 1]):
        fe_model.add_element(
            element_id=element_id + 1,
            nodes=ect[element_id, 2:].tolist(),
            part_id=int(part_id),
        )

    return fe_model


def map_to_part_id(
    centroid: np.ndarray, physical_space_map: np.ndarray, offset: int = 3
) -> int:
    """Map an element centroid to a part ID based on the spatial map.

    Args:
        centroid (np.ndarray): A 1D array representing the (x, y, z) coordinates of the element centroid.
        physical_space_map (np.ndarray): A 2D array where each row contains (x, y, z, part_id) for the spatial map.
        offset (int, optional): Offset to adjust part IDs. Defaults to 3.

    Returns:
        int: The part ID corresponding to the centroid.
    """
    # Create a KD tree to find the nearest neighbors
    physical_space_tree = sp.KDTree(physical_space_map[:, :3], leafsize=15)
    _, index = physical_space_tree.query(centroid, k=1)
    part_id = int(physical_space_map[index, 3]) + offset

    return part_id
