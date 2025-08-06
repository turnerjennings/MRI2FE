from typing import List, Optional

import numpy as np
import scipy.spatial as sp
from ants.core.ants_image import ANTsImage

from ..models.femodel import FEModel
from ..utilities import spatial_map


def map_MRE_to_mesh(
    mdl: FEModel,
    label_img: ANTsImage,
    region_properties: List,
    target_region_id: int = 4,
    label_background_id: int = 0,
    region_prefix: Optional[str] = None,
) -> FEModel:
    """Map MRE properties onto FE mesh and store associated material properties in the model material array

    Args:
        mdl (FEModel): FE model to map MRE regions onto
        label_img (ANTsImage): MRI image with integer labels for each MRE region
        region_properties (List): List of prony series properties for each MRE region
        target_region_id (int, optional): Target PID in the FE model to replace with segmented MRE IDs   .
        label_background_id (int, optional): Integer label in the label_img associated with the background. Defaults to 0.
        region_prefix (str, optional): name prefix for the new part assignment names

    Raises:
        TypeError: Input not the right type
        ValueError: negative target region

    Returns:
        FEModel: Model with updated PIDs for the target region and material properties added
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
        assert mdl.centroid_table is not None

    label_img_long = spatial_map(label_img)

    # find max ID in existing model to use as offset
    max_id = np.max(mdl.element_table[:, 1])

    # create filtered space map and centroids for ROI only
    label_img_long_nobackground = label_img_long[
        label_img_long[:, 3] != label_background_id
    ]

    # find elements and centroids within target label region
    ect_region_mask = mdl.element_table[:, 1] == target_region_id

    elcentroids_ROI = mdl.centroid_table[ect_region_mask, :]


    # create KDTree
    physical_space_tree = sp.KDTree(
        label_img_long_nobackground[:, 0:3], leafsize=15
    )

    d, idx = physical_space_tree.query(elcentroids_ROI, k=1)
    # calculate correct PID offset
    if max_id == target_region_id:
        offset = max_id - 1
    else:
        offset = max_id

    new_pids = label_img_long_nobackground[idx, 3] + offset

    mdl.element_table[ect_region_mask, 1] = new_pids

    # map part IDs and material ids

    for idx, mat in enumerate(region_properties):
        if idx is not label_background_id:
            mat_id = idx + offset

            # create part
            if region_prefix is not None:
                mdl.part_info[str(mat_id)] = {
                    "name": f"{region_prefix}_{idx}",
                    "constants": [mat_id],
                }
            else:
                mdl.part_info[str(mat_id)] = {
                    "name": f"{target_region_id}_{idx}",
                    "constants": [mat_id],
                }

            # create material
            if region_prefix is not None:
                mdl.material_info.append(
                    {
                        "type": "KELVIN-MAXWELL_VISCOELASTIC",
                        "name": f"{target_region_id}_{idx}",
                        "ID": mat_id,
                        "constants": [0, 0] + list(mat),
                    }
                )
            else:
                mdl.material_info.append(
                    {
                        "type": "KELVIN-MAXWELL_VISCOELASTIC",
                        "ID": mat_id,
                        "constants": [0, 0] + list(mat),
                    }
                )

    return mdl
