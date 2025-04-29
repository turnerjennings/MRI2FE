import ants
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
    DR_ants = ants.image_read(DR)

    SS_ants = ants.image_read(SS)

    bw_ants = ants.image_read(geom)
    bw_ants = ants.resample_image(bw_ants, SS_ants.shape, True)

    # segment MRI model using SS
    print(f"{datetime.now()}\t\tSegmenting models...")
    MRE_ants_segmented = ants.kmeans.kmeans_segmentation(SS_ants, n_segs)

    # calculate average values for SS and DR in each region
    print(f"{datetime.now()}\t\tCalculating segment means...")
    ROI_means = {"index": [], "Ginf": [], "G1": [], "Tau": []}

    for i in range(1, n_segs + 1):
        threshold = ants.threshold_image(
            MRE_ants_segmented["segmentation"], low_thresh=i, high_thresh=i + 1
        )

        SS_ROI = ants.mask_image(SS_ants, threshold)
        DR_ROI = ants.mask_image(DR_ants, threshold)

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

    tx = ants.registration(fixed=fi, moving=mi, type_of_transform="Elastic")

    # create plot outcome if requested
    if imgout is not None:
        ants.plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_sagittal.jpg",
            axis=0,
        )
        ants.plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_coronal.jpg",
            axis=1,
        )
        ants.plot(
            fi,
            overlay=mi,
            overlay_cmap="Dark2",
            overlay_alpha=0.8,
            filename=imgout + "_transverse.jpg",
            axis=2,
        )

    tx_output = tx["warpedmovout"]

    return tx_output, ROI_means


if __name__ == "__main__":
    BW_path = "ICBM/"
    BW_fname = "mni_icbm152_t1_tal_nlin_sym_09a.nii"

    mres = [
        5,
        37,
        49,
        51,
        61,
        67,
        78,
        91,
        93,
        124,
        10,
        17,
        20,
        28,
        66,
        75,
        81,
        89,
        99,
        128,
    ]
    ginf = []
    g1 = []
    tau = []
    mre_string = []

    for i in mres:
        MRE_path = "T:/LSDYNA/brainwebworkflow/MRE134"
        MRE_DR_fname = f"{i:03d}DR_warped.nii"
        MRE_SS_fname = f"{i:03d}Stiffness_warped.nii"

        imgout_path = "coregistration_verify"
        imgout_fname = (
            imgout_path + "/" + BW_fname[:-4] + "_" + MRE_DR_fname[0:3]
        )

        mapping, means = coregister_MRE(
            BW_path + "/" + BW_fname,
            MRE_path + "/" + MRE_DR_fname,
            MRE_path + "/" + MRE_SS_fname,
            imgout=imgout_fname,
        )

        # print(f"Ginf = {means['Ginf']}, G1={means['G1']}")
        ginf.append(means["Ginf"])
        g1.append(means["G1"])
        tau.append(means["Tau"])

    combined_data = np.column_stack((mres, ginf, g1, tau))
    header = "MRENo,Ginf1,Ginf2,Ginf3,Ginf4,Ginf5,G11,G12,G13,G14,G15,tau1,tau2,tau3,tau4,tau5"
    np.savetxt(
        "kmeans/MRE_stiffnesses.csv",
        combined_data,
        delimiter=",",
        header=header,
    )
