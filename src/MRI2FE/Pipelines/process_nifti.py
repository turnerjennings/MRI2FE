from lasso.dyna import D3plot
import numpy as np
import ants
from ..output.k_file_operations import parse_k_file, element_centroids
from ..Postprocess.calculate_strain import MRI_strain
from ..Postprocess.d3_to_nifti import (
    d3_to_displacement,
    save_field_variable,
    cloud_to_grid,
    grid_to_nifti,
)
from time import time
from typing import Literal
import datetime


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

            ux_img = ants.from_numpy(
                ux_grid,
                origin=template.origin,
                spacing=template.spacing,
                direction=template.direction,
            )

            uy = np.concatenate(
                (coords, disps[i, :, 1].reshape((disps.shape[1], 1))),
                axis=1,
            )

            uy_grid = cloud_to_grid(uy, dims=template.shape, lims=lims)

            uy_img = ants.from_numpy(
                uy_grid,
                origin=template.origin,
                spacing=template.spacing,
                direction=template.direction,
            )

            uz = np.concatenate(
                (coords, disps[i, :, 2].reshape((disps.shape[1], 1))),
                axis=1,
            )

            uz_grid = cloud_to_grid(uz, dims=template.shape, lims=lims)
            uz_img = ants.from_numpy(
                uz_grid,
                origin=template.origin,
                spacing=template.spacing,
                direction=template.direction,
            )

            strain_field = MRI_strain(
                ux=ux_img,
                uy=uy_img,
                uz=uz_img,
            )

            strain_field = ants.apply_transforms(
                fixed=icbm,
                moving=strain_field,
                transformlist=tx["fwdtransforms"],
                imagetype=1,
            )
            ants.image_write(
                strain_field, outpath + outname + f"_strainfield{i + 1}.nii"
            )
