from lasso.dyna import D3plot
import numpy as np
from scipy.interpolate import griddata
from ants.core.ants_image import ANTsImage
from ants import from_numpy, apply_transforms, image_write

from typing import Tuple


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


def grid_to_nifti(datagrid: np.ndarray, template: ANTsImage) -> ANTsImage:
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

    new_plot = from_numpy(
        data=datagrid, origin=origin, spacing=spacing, direction=direction
    )

    return new_plot


def save_field_variable(
    coordinates: np.ndarray,
    field_variable: np.ndarray,
    template: ANTsImage,
    icbm: ANTsImage,
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

    fv_ants_warped = apply_transforms(
        fixed=icbm, moving=fv_nifti, transformlist=tx, imagetype=0
    )

    # ants.plot(icbm,overlay=fv_ants_warped,overlay_alpha=0.7,overlay_cmap='viridis',axis=0)

    if save_plot:
        image_write(fv_ants_warped, fname)

    return fv_ants_warped
