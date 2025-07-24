import numpy as np
from ants.core.ants_image import ANTsImage
import scipy.spatial as sp
from ..utilities import spatial_map
from ..models.femodel import FEModel


def map_MRE_to_mesh(
    mdl: FEModel, label_img: ANTsImage, target_region_id: int = 4
) -> FEModel:
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
    if not isinstance(mdl, FEModel):
        raise TypeError("mdl must be a FEModel object")

    if not isinstance(label_img, ANTsImage):
        raise TypeError("map must be an ANTsImage object")

    if not isinstance(target_region_id, int):
        raise TypeError("offset must be an integer")
    if target_region_id < 0:
        raise ValueError("offset must be non-negative")

    # check for centroids
    if mdl.centroid_table is None:
        mdl.update_centroids()

    label_img_long = spatial_map(label_img)

    # find max ID in existing model to use as offset
    max_id = np.max(mdl.element_table[:, 1])

    # create filtered space map and centroids for ROI only
    label_img_long_nonzero = label_img_long[label_img_long[:, 3] > 0]

    # find elements and centroids within target label region
    ect_region_mask = mdl.element_table[:, 1] == target_region_id
    elcentroids_ROI = mdl.centroid_table[ect_region_mask, :]

    print(f"elements in the ROI: {elcentroids_ROI.shape}")
    print(f"points from the MRE label image: {label_img_long_nonzero.shape}")
    print(
        f"min/max of ROI: {np.min(elcentroids_ROI, axis=0)},{np.max(elcentroids_ROI, axis=0)}"
    )
    print(
        f"min/max of points in the MRE label image: {np.min(label_img_long_nonzero, axis=0)},{np.max(label_img_long_nonzero, axis=0)}"
    )

    # create KDTree
    physical_space_tree = sp.KDTree(
        label_img_long_nonzero[:, 0:3], leafsize=15
    )

    d, idx = physical_space_tree.query(elcentroids_ROI, k=1)
    # calculate correct PID offset
    if max_id == target_region_id:
        offset = max_id - 1
    else:
        offset = max_id

    print(f"idx shape: {idx.shape}\nunique indices found:{np.unique(idx)}")

    new_pids = label_img_long_nonzero[idx, 3] + offset
    print(f"min/max of new pids: {np.min(new_pids)},{np.max(new_pids)}")

    mdl.element_table[ect_region_mask, 1] = new_pids

    return mdl
