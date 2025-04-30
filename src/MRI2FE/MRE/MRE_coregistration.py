from ants import (
    image_read,
    resample_image,
    threshold_image,
    mask_image,
    registration,
    plot,
)
from ants import kmeans_segmentation
import numpy as np
from .calculate_prony import calculate_prony

from datetime import datetime


def coregister_MRE(
    geom: str, DR: str, SS: str, n_segs: int = 5, imgout: str = None
):
    """Generate MRE map segmented and transformed into the basis of a subject-specific geometry

    Args:
        geom (string): Filepath to the geometry MRI file
        DR (string): filepath to the damping ratio MRE file
        SS (string): filepath to the shear stiffness MRE file
        n_segs (int, optional): Number of segments to subdivide the MRE data. Defaults to 5.
        imgout (string, optional): filepath to save output image of final coregistration for verification

    Returns:
        tx_output (NiFTI1Image): MRE data co-registered onto geometry and segmented
        ROI_means (dict): table of prony series values for each index in the segmented model
    """

    # load MRI data
    print(f"{datetime.now()}\t\tLoading MRI...")

    # convert to ants format and resample
    print(f"{datetime.now()}\t\tConverting to ANTS format...")
    DR_ants = image_read(DR)

    SS_ants = image_read(SS)

    bw_ants = image_read(geom)
    bw_ants = resample_image(bw_ants, SS_ants.shape, True)

    # segment MRI model using SS
    print(f"{datetime.now()}\t\tSegmenting models...")
    MRE_ants_segmented = kmeans_segmentation(SS_ants, n_segs)

    # calculate average values for SS and DR in each region
    print(f"{datetime.now()}\t\tCalculating segment means...")
    ROI_means = {"index": [], "Ginf": [], "G1": [], "Tau": []}

    for i in range(1, n_segs + 1):
        threshold = threshold_image(
            MRE_ants_segmented["segmentation"], low_thresh=i, high_thresh=i + 1
        )

        SS_ROI = mask_image(SS_ants, threshold)
        DR_ROI = mask_image(DR_ants, threshold)

        SS_roi_np = SS_ROI.numpy()
        SS_roi_np_nonzero = SS_roi_np[SS_roi_np > 0.0]

        DR_roi_np = DR_ROI.numpy()
        DR_roi_np_nonzero = DR_roi_np[DR_roi_np > 0.0]

        SS_mean = np.mean(SS_roi_np_nonzero) / 1000
        DR_mean = np.mean(DR_roi_np_nonzero)

        Ginf, G1, tau, _, _ = calculate_prony(SS_mean, DR_mean, 50.0)

        print(
            f"Segment {i}: Mean SS = {SS_mean}, Mean DR = {DR_mean}, Ginf = {Ginf}, G1 = {G1}, tau = {tau}"
        )

        ROI_means["index"].append(i)
        ROI_means["Ginf"].append(Ginf / 1000)
        ROI_means["G1"].append(G1 / 1000)
        ROI_means["Tau"].append(tau)

    print(ROI_means)

    # coregister MRE onto BW geometry

    fi = bw_ants

    mi = MRE_ants_segmented["segmentation"]

    tx = registration(fixed=fi, moving=mi, type_of_transform="Elastic")

    # create plot outcome if requested
    if imgout is not None:
        plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_sagittal.jpg",
            axis=0,
        )
        plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_coronal.jpg",
            axis=1,
        )
        plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_transverse.jpg",
            axis=2,
        )

    tx_output = tx["warpedmovout"]

    return tx_output, ROI_means
