import numpy as np
import ants
import scipy.spatial as sp
from .k_file_operations import parse_k_file, element_centroids
from .MRE_coregistration import coregister_MRE
from .utilities import COM_align, ants_affine

import datetime


def spatial_map(infile: ants.core.ants_image.ANTsImage):
    """Convert data from NIFTI file to (4,n) array of x,y,z,voxel value in physical space

    Args:
        infile (nb.nifti1): NIFTI file to convert

    Returns:
        coordinates_and_values (np.ndarray): array of x,y,z coordinates in physical space with the corresponding voxel value
    """
    # extract data from NIFTI
    image_data = infile.numpy()
    affine = ants_affine(infile)
    dimensions = image_data.shape

    # create coordinates in voxel space
    x = np.arange(dimensions[0])
    y = np.arange(dimensions[1])
    z = np.arange(dimensions[2])
    xv, yv, zv = np.meshgrid(x, y, z, indexing="ij")

    voxel_coords = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
    voxel_coords_homogeneous = np.hstack(
        [voxel_coords, np.ones((voxel_coords.shape[0], 1))]
    )

    # map voxel coordinates to physical space
    physical_coords = voxel_coords_homogeneous @ affine.T
    voxel_values = np.round(image_data.ravel())
    coordinates_and_values = np.hstack(
        [physical_coords[:, :3], voxel_values[:, np.newaxis]]
    )

    return coordinates_and_values


def map_MRE_to_mesh(
    map: ants.core.ants_image.ANTsImage,
    elcentroids: np.ndarray,
    ect: np.ndarray,
    offset: int = 3,
    csvpath: str = None,
):
    """Map elements to parts from a segmented spatial map

    Args:
        map (AntsImage): segmented spatial map in voxel space
        elcentroids (np.ndarray): (3,n_elements) array of the coordinates of each element centroid
        ect (np.ndarray): element connectivity table in 10-node LS-Dyna format
        offset (int, optional): optional filter if it is not desired to remap all pids, will skip mapping any elements belonging to pid <= offset. Defaults to 3.

    Returns:
        ect (np.ndarray): new ect in 10-node LS-DYNA format with updated PIDs
    """

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

    return ect


if __name__ == "__main__":
    # load and map MRI data
    BW_path = "T:/LSDYNA/brainwebworkflow/brainwebmri"
    BW_fname = "subject04_crisp_v.mnc"

    MRE_path = "T:/LSDYNA/brainwebworkflow/MRE134"
    MRE_DR_fname = "001DR_warped.nii"
    MRE_SS_fname = "001Stiffness_warped.nii"

    mapping, means = coregister_MRE(
        BW_path + "/" + BW_fname,
        MRE_path + "/" + MRE_DR_fname,
        MRE_path + "/" + MRE_SS_fname,
    )

    # load and process k file
    fepath = "brainwebFE/BrainWeb_Subject04_updated.k"
    print("processing k file...")
    ect_array, node_array = parse_k_file(fepath)

    print(
        f"FE model: min/max x=({np.min(node_array[:, 1])},{np.max(node_array[:, 1])}), y=({np.min(node_array[:, 2])},{np.max(node_array[:, 2])}), z=({np.min(node_array[:, 3])},{np.max(node_array[:, 3])})"
    )
    print("calculating centroids...")
    centroids = np.apply_along_axis(
        element_centroids, 1, ect_array, node_array
    )

    print("Mapping MRE to mesh...")
    ect_array_mapped = map_MRE_to_mesh(mapping, centroids, ect_array)
