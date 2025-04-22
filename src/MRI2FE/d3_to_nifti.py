from lasso.dyna import D3plot
import numpy as np
from scipy.interpolate import griddata
import ants
from .k_file_operations import parse_k_file, element_centroids
from .utilities import get_spacing
from time import time
from typing import Literal, Tuple

import datetime


def calculate_strain(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    spacing: Tuple,
    stype: Literal["eng", "green"] = "green",
    print_summary: bool = False,
) -> np.ndarray:
    """Function to calculate a 3d strain field given x,y,z displacements"""

    if not (ux.shape == uy.shape == uz.shape):
        raise IndexError("ux,uy,uz not the same shape")

    gradients = [np.gradient(u, *spacing) for u in [ux, uy, uz]]

    # Initialize the strain tensor array
    shape = ux.shape[:3] + (3, 3)
    strain_tensor = np.zeros(shape)

    if stype == "green":
        # BUG massive overcalculation here
        # Compute the strain tensor
        for i in range(3):  # Loop over displacement components (ux, uy, uz)
            for j in range(3):  # Loop over spatial dimensions (x, y, z)
                strain_tensor[..., i, j] = 0.5 * (
                    gradients[i][j]
                    + gradients[j][i]
                    + gradients[0][i] * gradients[0][i]
                    + gradients[1][i] * gradients[1][i]
                    + gradients[2][i] * gradients[2][i]
                )
    elif stype == "eng":
        # Compute the strain tensor
        for i in range(3):  # Loop over displacement components (ux, uy, uz)
            for j in range(3):  # Loop over spatial dimensions (x, y, z)
                strain_tensor[..., i, j] = 0.5 * (
                    gradients[i][j] + gradients[j][i]
                )
    else:
        raise ValueError("stype must be 'green' or 'eng'")

    strain_stack = np.stack(
        [
            strain_tensor[:, :, :, 0, 0],
            strain_tensor[:, :, :, 1, 1],
            strain_tensor[:, :, :, 2, 2],
            strain_tensor[:, :, :, 0, 1],
            strain_tensor[:, :, :, 1, 2],
            strain_tensor[:, :, :, 0, 2],
        ],
        axis=-1,
    )

    strain_stack = np.nan_to_num(strain_stack)
    return strain_stack


def d3_to_displacement(
    plot: D3plot, part_filter: list = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Take d3eigv file as input and output the node coordinates and displacements for relevant parameters

    Args:
        plot (D3plot): D3eigv file for model
        part_filter (list, optional): Subset of parts to evaluate displacements at. Defaults to None.

    Returns:
        coords_filtered (np.ndarray): node coordinates (reference config), (n_nodes, 3)
        disp_filtered (np.ndarray): node displacements (deformed config - reference config), (n_timesteps, n_nodes, 3)
    """

    # find nodes associated with part_filter
    node_indexes = plot.arrays["element_solid_node_indexes"]
    if part_filter is not None:
        part_ids = plot.arrays["element_solid_part_indexes"]

        element_mask = np.isin(part_ids, part_filter)

        nodes_to_eval = np.unique(node_indexes[element_mask, :])

    else:
        nodes_to_eval = np.unique(node_indexes)

    # print(nodes_to_eval.shape)

    # extract node displacement
    disp = plot.arrays["node_displacement"]
    # print(f"disp shape: {disp.shape}")
    coords = plot.arrays["node_coordinates"]

    # filter coords and displacements
    coords_filtered = coords[nodes_to_eval, :]
    # coords_filtered[:,[1,2]] = coords_filtered[:,[2,1]]

    disp_filtered = disp[:, nodes_to_eval, :]
    # disp_filtered[:,[1,2]] = disp[:,[2,1]]

    return coords_filtered, disp_filtered - coords_filtered


def cloud_to_grid(
    data: np.ndarray, dims: tuple, lims: dict = None
) -> np.ndarray:
    """Convert point cloud data (Lagrangian coordinates) to structured grid data (Eulerian coordinates)

    Args:
        Data (np.ndarray): point cloud coordinates and point values in format (x,y,z,val)
        dims (tuple[int]): voxel cloud dimensions for grid coordinate field
        lims (dict, optional): optional min/max spatial coordinate limits for voxel field.
            calculated from min/max coordinates otherwise.

    Returns:
        data_out (np.ndarray): grid with shape (dims) linearly interpolated with point cloud data
    """

    # find minima and maxima
    if lims is None:
        minima = np.min(data, axis=0)
        maxima = np.max(data, axis=0)
    else:
        minima = lims["min"]
        maxima = lims["max"]

    # print(f"minima: {minima}, maxima: {maxima}")

    # create grid
    x = np.linspace(minima[0], maxima[0], dims[0])
    y = np.linspace(minima[1], maxima[1], dims[1])
    z = np.linspace(minima[2], maxima[2], dims[2])

    # print(f"xmin: {np.min(x)},xmax: {np.max(x)}")

    xv, yv, zv = np.meshgrid(x, y, z, indexing="ij")
    # print(f"xvmin: {np.min(xv)},xvmax: {np.max(xv)}")
    grid_coords = np.column_stack((xv.ravel(), yv.ravel(), zv.ravel()))
    # np.savetxt('test_output/grid_coords.csv',grid_coords,delimiter=',',header="x,y,z,xdisp")

    # print(f"grid coords shape: {grid_coords.shape}")

    # map data onto grid
    data_out = griddata(
        data[:, 0:3], data[:, 3], grid_coords, method="linear", fill_value=0.0
    )
    # print(f"grid coords shape: {grid_coords.shape}, data out shape:{data_out.shape}")
    # write_output = np.column_stack((grid_coords,data_out))
    # np.savetxt('test_output/mode_1_xdisp_grid.csv',write_output,delimiter=',',header="x,y,z,xdisp")
    data_out = data_out.reshape(dims)

    return data_out


def grid_to_nifti(
    datagrid: np.ndarray, template: ants.core.ants_image.ANTsImage
) -> ants.core.ants_image.ANTsImage:
    """Convert an input numpy array voxel grid to a nifti file based on a template input nifti

    Args:
        datagrid (np.ndarray): Voxel grid of data
        template (np.ndarray): template nifti file with header and affine transformation info

    Returns:
        new_plot (AntsImage): nifti file containing voxel data
    """

    spacing = template.spacing
    direction = template.direction
    origin = template.origin

    new_plot = ants.from_numpy(
        data=datagrid, origin=origin, spacing=spacing, direction=direction
    )

    return new_plot


def save_field_variable(
    coordinates: np.ndarray,
    field_variable: np.ndarray,
    template: ants.core.ants_image.ANTsImage,
    icbm: ants.core.ants_image.ANTsImage,
    tx: str,
    fname: str,
    step: int = 0,
    limits: dict = {"max": [1, 1, 1], "min": [0, 0, 0]},
    shape: Tuple = None,
    save_plot: bool = True,
) -> None:
    """Function to map a field variable to a 3d grid, transform into ICBM space, and save to file.

    Args:
        coordinates (np.ndarray): Array of shape (n_points, 3) with x,y,z coordinates for each point to evaluate.
        field_variable (np.ndarray): Array of shape (n_timesteps,n_points) with field variable values.
        template (AntsImage): Nifti image with spacing, direction, transform information
        icbm (AntsImage): Nifti image with ICBM atlas model
        tx (str): Ants transform mapping the original data to the ICBM model
        fname (str): File path to save result to
        step (int): Timestep to save field variable from
        limits (dict): dict with two keys (min and max) storing the spatial limits of the model
        shape (tuple): tuple with desired spatial resolution for the model.
    """

    if coordinates.shape[0] != field_variable.shape[1]:
        raise IndexError(
            "Number of points on axis 0 of coordinates does not equal number of points on axis 1 of field_variable"
        )

    print(
        f"coordinates shape: {coordinates.shape}, field variable shape: {field_variable.shape}"
    )

    if shape is None:
        shape = template.get_fdata().shape

    coord_fv_stack = np.column_stack((coordinates, field_variable[step, :]))

    # convert stacked values to ants image
    fv_grid = cloud_to_grid(coord_fv_stack, shape, limits)
    fv_nifti = grid_to_nifti(fv_grid, template)

    # check if scalar field or 4d field

    fv_ants_warped = ants.apply_transforms(
        fixed=icbm, moving=fv_nifti, transformlist=tx, imagetype=0
    )

    # ants.plot(icbm,overlay=fv_ants_warped,overlay_alpha=0.7,overlay_cmap='viridis',axis=0)

    if save_plot:
        ants.image_write(fv_ants_warped, fname)

    return fv_ants_warped


def process_nifti(
    icbm: ants.core.ants_image.ANTsImage,
    d3: D3plot,
    mdl: str,
    template: ants.core.ants_image.ANTsImage,
    outpath: str,
    outname: str,
    output_type: Literal["displacement", "strain"] = "displacement",
    dims=(False, False, False, True),
    save_states: bool = True,
    states_to_save: Literal["all", "max"] = "all",
):
    """Save displacement data from a d3 state file to a series of nifti files co-registered to the ICBM-152 atlas.

    Args:
        icbm (ANTsImage): Nifti structure containing ICBM-152 T1-weighted atlas
        d3 (D3plot): d3plot/eigv/ssd/etc structure containing the data to save
        mdl (str): file path to the original model geometry keyword file
        template (ANTsImage): Template nifti file with header and affine transformation info
        outpath (str): File path to save results to
        outname (str): File name root, to be appended with data direction
        dims (tuple, optional): Which data dimensions to save, (x,y,z,magnitude). Defaults to (False, False, False, True).
        save_states (bool, optional): Whether to save resultant plots to files, defaults to True
        states_to_save ("max" or "all", optional): If "all" selected, all states in the d3 file will be saved.  If "max" selected, only the state with highest average deformation magnitude will be saved. Defaults to "all".
    """

    # calculate spatial limits of model
    step_start = time()
    print(f"{datetime.now()}\tCalculating model limits...")
    # load coords and fix directions
    d3.arrays["node_coordinates"][:, [1, 0]] = d3.arrays["node_coordinates"][
        :, [0, 1]
    ]
    d3.arrays["node_displacement"][:, :, [1, 0]] = d3.arrays[
        "node_displacement"
    ][:, :, [0, 1]]

    lims = {
        "min": np.min(d3.arrays["node_coordinates"], axis=0),
        "max": np.max(d3.arrays["node_coordinates"], axis=0),
    }

    step_end = time()
    print(
        f"{datetime.now()}\tComplete... {step_end - step_start:.2f} seconds."
    )

    freqs = d3.arrays["timesteps"]

    # create nifti file with pids
    step_start = time()
    print(f"{datetime.now()}\tCreating nifti for part IDs...")
    ect_array, node_array = parse_k_file(mdl)

    centroids = np.apply_along_axis(
        element_centroids, 1, ect_array, node_array
    )
    centroids[:, [0, 1]] = centroids[:, [1, 0]]

    pids = d3.arrays["element_solid_part_indexes"]

    pids = np.column_stack((centroids, pids))

    pid_grid = cloud_to_grid(pids, template.shape, lims)

    pid_plot = grid_to_nifti(pid_grid, template)

    # ants.plot(icbm, overlay=pid_plot)
    print("Saving nifti for PIDs...")
    step_end = time()
    print(
        f"{datetime.now()}\tComplete... {step_end - step_start:.2f} seconds."
    )

    step_start = time()
    print(f"{datetime.now()}\tCoregistering part IDs to ICBM template...")

    # ants.plot(icbm,overlay=pid_plot)

    mi = ants.threshold_image(pid_plot, low_thresh=3.0)
    fi = ants.threshold_image(icbm, low_thresh=65, high_thresh=100)
    tx = ants.registration(fixed=fi, moving=mi, type_of_transform="SyN")

    # ants.plot(icbm,overlay=tx["warpedmovout"])
    ants.image_write(tx["warpedmovout"], outpath + outname + "_pids.nii")

    transform = tx["fwdtransforms"]

    step_end = time()
    print(
        f"{datetime.now()}\tComplete... {step_end - step_start:.2f} seconds."
    )

    print(f"{datetime.now()}\tGenerating displacement grids...")
    coords, disps = d3_to_displacement(d3, [3, 4, 5, 6, 7])
    # print(f"coords shape: {coords.shape}, disps shape: {disps.shape}")
    disps_mag = np.linalg.norm(disps, axis=2)
    # print(f"disps_mag shape: {disps_mag.shape}")
    if states_to_save == "all":
        states = range(freqs.shape[0])

    elif states_to_save == "max":
        disps_mag_avg = np.mean(disps_mag, axis=1)
        max_idx = np.argmax(disps_mag_avg)
        states = [max_idx]

    # create x,y,z,mag nifti for each mode
    for i in states:
        # check if output is displacement
        if freqs[i] > 1 and output_type == "displacement":
            print(f"{datetime.now()}\tProcessing mode {i + 1}...")
            step_start = time()

            if dims[0]:
                print(
                    f"{datetime.now()}\t\tmapping x-disp to nifti for mode {i + 1}..."
                )
                _ = save_field_variable(
                    coordinates=coords,
                    field_variable=disps[:, :, 0],
                    template=template,
                    icbm=icbm,
                    tx=transform,
                    fname=outpath + outname + f"_mode{i + 1}x.nii",
                    step=i,
                    limits=lims,
                    shape=template.shape,
                )

            if dims[1]:
                print(
                    f"{datetime.now()}\t\tmapping y-disp to nifti for mode {i + 1}..."
                )
                _ = save_field_variable(
                    coordinates=coords,
                    field_variable=disps[:, :, 1],
                    template=template,
                    icbm=icbm,
                    tx=transform,
                    fname=outpath + outname + f"_mode{i + 1}y.nii",
                    step=i,
                    limits=lims,
                    shape=template.shape,
                )

            if dims[2]:
                print(
                    f"{datetime.now()}\t\tmapping z-disp to nifti for mode {i + 1}..."
                )
                _ = save_field_variable(
                    coordinates=coords,
                    field_variable=disps[:, :, 2],
                    template=template,
                    icbm=icbm,
                    tx=transform,
                    fname=outpath + outname + f"_mode{i + 1}z.nii",
                    step=i,
                    limits=lims,
                    shape=template.shape,
                )

            if dims[3]:
                print(
                    f"{datetime.now()}\t\tmapping mag-disp to nifti for mode {i + 1}..."
                )
                _ = save_field_variable(
                    coordinates=coords,
                    field_variable=disps_mag,
                    template=template,
                    icbm=icbm,
                    tx=transform,
                    fname=outpath + outname + f"_mode{i + 1}mag.nii",
                    step=i,
                    limits=lims,
                    shape=template.shape,
                )

            step_end = time()
            print(
                f"{datetime.now()}\tProcessing mode {i + 1} complete... {step_end - step_start:.2f} seconds."
            )

        # strain output type
        elif freqs[i] > 0 and output_type == "strain":
            print(f"{datetime.now()}\tProcessing strain for step {i + 1}...")
            step_start = time()

            ux = np.concatenate(
                (coords, disps[i, :, 0].reshape((disps.shape[1], 1))),
                axis=1,
            )

            ux_grid = cloud_to_grid(ux, dims=template.shape, lims=lims)

            uy = np.concatenate(
                (coords, disps[i, :, 1].reshape((disps.shape[1], 1))),
                axis=1,
            )

            uy_grid = cloud_to_grid(uy, dims=template.shape, lims=lims)

            uz = np.concatenate(
                (coords, disps[i, :, 2].reshape((disps.shape[1], 1))),
                axis=1,
            )

            uz_grid = cloud_to_grid(uz, dims=template.shape, lims=lims)

            strain_field = calculate_strain(
                ux=ux_grid,
                uy=uy_grid,
                uz=uz_grid,
                spacing=get_spacing(lims=lims, dims=template.shape),
                stype="eng",
                print_summary=True,
            )

            strain_ants = ants.from_numpy(
                data=strain_field,
                origin=template.origin,
                spacing=template.spacing,
                direction=template.direction,
                has_components=True,
            )

            strain_ants = ants.apply_transforms(
                fixed=icbm,
                moving=strain_ants,
                transformlist=tx["fwdtransforms"],
                imagetype=1,
            )
            ants.image_write(
                strain_ants, outpath + outname + f"_strainfield{i + 1}.nii"
            )


if __name__ == "__main__":
    niftipath = "t:\\LSDYNA\\brainwebworkflow\\transforms\\0.nii"
    plotpath = "t:\\LSDYNA\\brainwebworkflow\\test_output\\d3eigv"
