from ants import (
    image_read,
    resample_image,
    threshold_image,
    mask_image,
    registration,
    plot,
    kmeans_segmentation,
)
import numpy as np
from .calculate_prony import calculate_prony
from datetime import datetime

def coregister_MRE_images(geom: str, DR: str, SS: str, imgout: str = None):
    """Coregister MRE images to MRI geometry.

    Args:
        geom (str): Filepath to the geometry MRI file.
        DR (str): Filepath to the damping ratio MRE file.
        SS (str): Filepath to the shear stiffness MRE file.
        imgout (str, optional): Filepath prefix to save visualization images.

    Returns:
        dict: Contains coregistered SS, DR images and transformation metadata.
    """
    print(f"{datetime.now()}\t\tLoading MRI and MRE...")
    DR_ants = image_read(DR)
    SS_ants = image_read(SS)
    geom_ants = image_read(geom)
    geom_ants = resample_image(geom_ants, SS_ants.shape, True)

    print(f"{datetime.now()}\t\tPerforming registration...")
    tx = registration(fixed=geom_ants, moving=SS_ants, type_of_transform="Elastic")

    SS_registered = tx["warpedmovout"]
    DR_registered = registration(fixed=geom_ants, moving=DR_ants,
                                  type_of_transform="Elastic")["warpedmovout"]

    if imgout is not None:
        plot(geom_ants, overlay=SS_registered, overlay_cmap="Dark2", overlay_alpha=0.8,
             filename=imgout + "_sagittal.jpg", axis=0)
        plot(geom_ants, overlay=SS_registered, overlay_cmap="Dark2", overlay_alpha=0.8,
             filename=imgout + "_coronal.jpg", axis=1)
        plot(geom_ants, overlay=SS_registered, overlay_cmap="Dark2", overlay_alpha=0.8,
             filename=imgout + "_transverse.jpg", axis=2)

    return {"SS": SS_registered, "DR": DR_registered, "tx": tx}


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
        print(f"Segment {i}: SS={SS_mean:.3f}, DR={DR_mean:.3f}, Ginf={Ginf:.3f}, G1={G1:.3f}, tau={tau:.3f}")

        ROI_means["index"].append(i)
        ROI_means["Ginf"].append(Ginf / 1000)
        ROI_means["G1"].append(G1 / 1000)
        ROI_means["Tau"].append(tau)

    return ROI_means

def run_MRE_pipeline(geom, DR, SS, n_segs=5, imgout=None):
    coreg = coregister_MRE_images(geom, DR, SS, imgout)
    ROI_means = segment_MRE_regions(coreg["SS"], coreg["DR"], n_segs)
    return coreg["SS"], ROI_means