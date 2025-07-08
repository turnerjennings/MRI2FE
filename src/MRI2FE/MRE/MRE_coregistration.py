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
from typing import Union, List
import os

def coregister_MRE_images(
    geom: Union[str, ANTsImage],
    geom_mask: Union[str, ANTsImage] = None,
    gp_list: List[Union[str, ANTsImage]] = None,
    gpp_list: List[Union[str, ANTsImage]] = None,
    mu_list: List[Union[str, ANTsImage]] = None,
    xi_list: List[Union[str, ANTsImage]] = None,
    imgout: str = None,
):
    """Coregister multiple MRE images to MRI geometry.

    Args:
        geom: Image or filepath to the geometry MRI file.
        geom_mask: Image or filepath to the geometry mask.
        gp_list: List of storage modulus maps.
        gpp_list: List of loss modulus maps.
        mu_list: List of mu parameter maps.
        xi_list: List of xi parameter maps.
        imgout: Output prefix for saved visualizations.

    Returns:
        List of dicts with coregistered results and transformation metadata.
    """
    if geom is None:
        raise ValueError("Geometry image is required")

    # Load geometry image
    if isinstance(geom, str):
        if not os.path.exists(geom):
            raise FileNotFoundError(f"Geometry image file not found: {geom}")
        geom = image_read(geom)
    elif not isinstance(geom, ANTsImage):
        raise TypeError("geom must be either a filepath string or ANTsImage object")

    # Load optional geometry mask
    if geom_mask is not None:
        if isinstance(geom_mask, str):
            if not os.path.exists(geom_mask):
                raise FileNotFoundError(f"Geometry mask file not found: {geom_mask}")
            geom_mask = image_read(geom_mask)
        elif not isinstance(geom_mask, ANTsImage):
            raise TypeError("geom_mask must be either a filepath string or ANTsImage object")

    results = []

    # Coregistration for complex shear modulus (gp + gpp)
    if gp_list and gpp_list:
        for idx, (gp, gpp) in enumerate(zip(gp_list, gpp_list)):
            # Load images
            if isinstance(gp, str):
                gp = image_read(gp)
            if isinstance(gpp, str):
                gpp = image_read(gpp)

            # Create binary mask for moving image
            gp_mask = threshold_image(gp, np.min(np.nonzero(gp.numpy())), np.max(np.nonzero(gp.numpy())))

            # Resample geometry to match moving image resolution
            geom_ants = resample_image_to_target(geom, gp)

            # Perform registration
            tx = registration(
                fixed=geom_ants,
                moving=gp,
                moving_mask=gp_mask,
                mask=geom_mask,
                type_of_transform="Elastic",
            )

            # Apply transformation to both gp and gpp
            gp_out = tx["warpedmovout"]
            gpp_out = apply_transforms(fixed=geom, moving=gpp, transformlist=tx["fwdtransforms"])

            # Save plots if requested
            if imgout is not None:
                base = f"{imgout}_pair{idx}"
                plot(geom, overlay=gp_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_sagittal.jpg", axis=0)
                plot(geom, overlay=gp_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_coronal.jpg", axis=1)
                plot(geom, overlay=gp_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_transverse.jpg", axis=2)

            results.append({"gp": gp_out, "gpp": gpp_out, "transform": tx["fwdtransforms"]})

    # Coregistration for shear stiffness and damping ratio (mu + xi)
    elif mu_list and xi_list:
        for idx, (mu, xi) in enumerate(zip(mu_list, xi_list)):
            # Load images
            if isinstance(mu, str):
                mu = image_read(mu)
            if isinstance(xi, str):
                xi = image_read(xi)

            # Create binary mask for moving image
            mu_mask = threshold_image(mu, np.min(np.nonzero(mu.numpy())), np.max(np.nonzero(mu.numpy())))

            # Resample geometry to match moving image resolution
            geom_ants = resample_image_to_target(geom, mu)

            # Perform registration
            tx = registration(
                fixed=geom_ants,
                moving=mu,
                moving_mask=mu_mask,
                mask=geom_mask,
                type_of_transform="Elastic",
            )

            # Apply transformation to both mu and xi
            mu_out = tx["warpedmovout"]
            xi_out = apply_transforms(fixed=geom, moving=xi, transformlist=tx["fwdtransforms"])

            # Save plots if requested
            if imgout is not None:
                base = f"{imgout}_pair{idx}"
                plot(geom, overlay=mu_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_sagittal.jpg", axis=0)
                plot(geom, overlay=mu_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_coronal.jpg", axis=1)
                plot(geom, overlay=mu_out, overlay_cmap="Dark2", overlay_alpha=0.8, filename=base + "_transverse.jpg", axis=2)

            results.append({"mu": mu_out, "xi": xi_out, "transform": tx["fwdtransforms"]})

    else:
        raise ValueError("Must provide either (gp_list, gpp_list) or (mu_list, xi_list)")

    return results

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
        # Threshold to create mask for each region
        threshold = threshold_image(segments, low_thresh=i, high_thresh=i + 1)
        SS_ROI = mask_image(SS_img, threshold)
        DR_ROI = mask_image(DR_img, threshold)

        # Extract non-zero voxels and calculate means
        SS_vals = SS_ROI.numpy()[SS_ROI.numpy() > 0.0]
        DR_vals = DR_ROI.numpy()[DR_ROI.numpy() > 0.0]

        SS_mean = np.mean(SS_vals) / 1000
        DR_mean = np.mean(DR_vals)

        # Compute Prony series parameters
        Ginf, G1, tau, _, _ = calculate_prony(SS_mean, DR_mean, 50.0)
        print(f"Segment {i}: SS={SS_mean:.3f}, DR={DR_mean:.3f}, Ginf={Ginf:.3f}, G1={G1:.3f}, tau={tau:.3f}")

        ROI_means["index"].append(i)
        ROI_means["Ginf"].append(Ginf / 1000)
        ROI_means["G1"].append(G1 / 1000)
        ROI_means["Tau"].append(tau)

    return ROI_means

def run_MRE_pipeline(geom, DR_list, SS_list, n_segs=5, imgout=None):
    """Run full MRE pipeline: coregistration + segmentation + Prony analysis.

    Args:
        geom: Geometry MRI file (str or ANTsImage).
        DR_list: List of damping ratio (G'') images.
        SS_list: List of shear stiffness (G') images.
        n_segs (int): Number of segments for clustering.
        imgout: Optional prefix for saving visualization images.

    Returns:
        Tuple of coregistration results and ROI statistics per image pair.
    """
    coreg_results = coregister_MRE_images(geom=geom, gp_list=SS_list, gpp_list=DR_list, imgout=imgout)
    all_results = []
    for result in coreg_results:
        ROI_means = segment_MRE_regions(result["gp"], result["gpp"], n_segs)
        all_results.append(ROI_means)
    return coreg_results, all_results
