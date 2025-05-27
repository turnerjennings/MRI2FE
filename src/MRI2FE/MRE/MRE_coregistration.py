from ants import (
    image_read,
    resample_image_to_target,
    threshold_image,
    mask_image,
    registration,
    plot,
    kmeans_segmentation,
    apply_transforms,
)
from ants.core.ants_image import ANTsImage
import numpy as np
from .calculate_prony import calculate_prony
from datetime import datetime
from typing import Union
import os


def coregister_MRE_images(
    geom: Union[str, ANTsImage],
    geom_mask: Union[str, ANTsImage] = None,
    gp: Union[str, ANTsImage] = None,
    gpp: Union[str, ANTsImage] = None,
    mu: Union[str, ANTsImage] = None,
    xi: Union[str, ANTsImage] = None,
    imgout: str = None,
):
    """Coregister MRE images to MRI geometry.

    Args:
        geom (str or ANTsImage): Image or filepath to the geometry MRI file.
        geom_mask (str or ANTsImage, optional): Image or filepath to the geometry mask.
        gp (str or ANTsImage): Image or filepath to the storage modulus map.
        gpp (str or ANTsImage): Image or filepath to the loss modulus map.
        mu (str or ANTsImage): Image or filepath to the mu parameter map.
        xi (str or ANTsImage): Image or filepath to the xi parameter map.
        imgout (str, optional): Filepath prefix to save visualization images.

    Raises:
        ValueError: incorrect input formats
        TypeError: If input types are invalid
        FileNotFoundError: If image files don't exist

    Returns:
        dict: Contains coregistered SS, DR images and transformation metadata.
    """
    if geom is None:
        raise ValueError("Geometry image is required")

    # Validate and load geometry image
    if isinstance(geom, str):
        if not os.path.exists(geom):
            raise FileNotFoundError(f"Geometry image file not found: {geom}")
        geom = image_read(geom)
    elif not isinstance(geom, ANTsImage):
        raise TypeError("geom must be either a filepath string or ANTsImage object")

    # Validate and load geometry mask if provided
    if geom_mask is not None:
        if isinstance(geom_mask, str):
            if not os.path.exists(geom_mask):
                raise FileNotFoundError(f"Geometry mask file not found: {geom_mask}")
            geom_mask = image_read(geom_mask)
        elif not isinstance(geom_mask, ANTsImage):
            raise TypeError("geom_mask must be either a filepath string or ANTsImage object")
        
        # Check mask dimensions match geometry
        if geom_mask.dimension != geom.dimension:
            raise ValueError(f"Geometry mask dimensions ({geom_mask.dimension}) do not match geometry image dimensions ({geom.dimension})")

    # Validate output path if provided
    if imgout is not None:
        output_dir = os.path.dirname(imgout)
        if output_dir and not os.path.exists(output_dir):
            raise ValueError(f"Output directory does not exist: {output_dir}")

    # check for valid input parameter sets
    first_set_valid = all(param is not None for param in [geom, gp, gpp])
    second_set_valid = all(param is not None for param in [geom, mu, xi])
    
    # Raise error if neither parameter set is complete
    if not (first_set_valid or second_set_valid):
        raise ValueError(
            "You must provide either (gp, gpp, w) or (mu, xi, w) as inputs"
        )

    # load geometry and mask (if applicable)
    if type(geom) is str:
        geom = image_read(geom)MRE

    if geom_mask is not None and type(geom_mask) is str:
        geom_mask = image_read(geom_mask)

    # if complex shear modulus is provided
    if first_set_valid:
        # load images
        if type(gp) is str:
            gp = image_read(gp)

        if type(gpp) is str:
            gpp = image_read(gpp)

        # find threshold limits and generate binary mask
        gp_arr = np.nonzero(gp.numpy())
        gp_min = np.min(gp_arr)
        gp_max = np.max(gp_arr)
        gp_mask = threshold_image(gp, gp_min, gp_max)

        # generate transformation
        geom_ants = resample_image_to_target(geom, gp)

        if geom_mask is not None:
            tx = registration(
                fixed=geom_ants,
                mask=geom_mask,
                moving=gp,
                moving_mask=gp_mask,
                type_of_transform="Elastic",
            )

        else:
            tx = registration(
                fixed=geom_ants,
                moving=gp,
                moving_mask=gp_mask,
                type_of_transform="Elastic",
            )

        gp_out = tx["warpedmovout"]

        transform = tx["fwdtransforms"]

        gpp_out = apply_transforms(
            fixed=geom, moving=gpp, transformlist=transform
        )

        out_dict = {"gp": gp_out, "gpp": gpp_out, "transform": transform}

        # write images to file if requested
        if imgout is not None:
            plot(
                geom,
                overlay=gp_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_sagittal.jpg",
                axis=0,
            )
            plot(
                geom,
                overlay=gp_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_coronal.jpg",
                axis=1,
            )
            plot(
                geom,
                overlay=gp_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_transverse.jpg",
                axis=2,
            )

        return out_dict

    # if ss/dr is provided
    else:
        # load images
        if type(mu) is str:
            mu = image_read(mu)

        if type(xi) is str:
            xi = image_read(xi)

        # find threshold limits and generate binary mask
        mu_arr = np.nonzero(mu.numpy())
        mu_min = np.min(mu_arr)
        mu_max = np.max(mu_arr)
        mu_mask = threshold_image(mu, mu_min, mu_max)

        # generate transformation
        geom_ants = resample_image_to_target(geom, mu)

        if geom_mask is not None:
            tx = registration(
                fixed=geom_ants,
                mask=geom_mask,
                moving=mu,
                moving_mask=mu_mask,
                type_of_transform="Elastic",
            )

        else:
            tx = registration(
                fixed=geom_ants,
                moving=mu,
                moving_mask=mu_mask,
                type_of_transform="Elastic",
            )

        mu_out = tx["warpedmovout"]

        transform = tx["fwdtransforms"]

        xi_out = apply_transforms(
            fixed=geom, moving=xi, transformlist=transform
        )

        out_dict = {"mu": mu_out, "xi": xi_out, "transform": transform}

        # write images to file if requested
        if imgout is not None:
            plot(
                geom,
                overlay=mu_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_sagittal.jpg",
                axis=0,
            )
            plot(
                geom,
                overlay=mu_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_coronal.jpg",
                axis=1,
            )
            plot(
                geom,
                overlay=mu_out,
                overlay_cmap="Dark2",
                overlay_alpha=0.8,
                filename=imgout + "_transverse.jpg",
                axis=2,
            )

        return out_dict


def segment_MRE_regions(SS_img, DR_img, n_segs: int = 5):
    """Segment MRE images and calculate Prony model parameters for each region.

    Args:
        SS_img: ANTsImage of shear stiffness.
        DR_img: ANTsImage of damping ratio.
        n_segs (int): Number of segments.

    Returns:
        dict: Dictionary of segment means and Prony parameters.
    """
    print(f"{datetime.now()}\t\tSegmenting MRE image...")
    segments = kmeans_segmentation(SS_img, n_segs)["segmentation"]
    ROI_means = {"index": [], "Ginf": [], "G1": [], "Tau": []}

    print(f"{datetime.now()}\t\tCalculating segment means...")
    for i in range(1, n_segs + 1):
        threshold = threshold_image(segments, low_thresh=i, high_thresh=i + 1)
        SS_ROI = mask_image(SS_img, threshold)
        DR_ROI = mask_image(DR_img, threshold)

        SS_vals = SS_ROI.numpy()[SS_ROI.numpy() > 0.0]
        DR_vals = DR_ROI.numpy()[DR_ROI.numpy() > 0.0]

        SS_mean = np.mean(SS_vals) / 1000
        DR_mean = np.mean(DR_vals)

        Ginf, G1, tau, _, _ = calculate_prony(SS_mean, DR_mean, 50.0)
        print(
            f"Segment {i}: SS={SS_mean:.3f}, DR={DR_mean:.3f}, Ginf={Ginf:.3f}, G1={G1:.3f}, tau={tau:.3f}"
        )

        ROI_means["index"].append(i)
        ROI_means["Ginf"].append(Ginf / 1000)
        ROI_means["G1"].append(G1 / 1000)
        ROI_means["Tau"].append(tau)

    return ROI_means


def run_MRE_pipeline(geom, DR, SS, n_segs=5, imgout=None):
    coreg = coregister_MRE_images(geom, DR, SS, imgout)
    ROI_means = segment_MRE_regions(coreg["SS"], coreg["DR"], n_segs)
    return coreg["SS"], ROI_means
