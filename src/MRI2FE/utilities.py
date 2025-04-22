import numpy as np
import ants
from datetime import datetime
from typing import Union
import re
import os


def create_title(bw_fname: str, MRE_fname: str):
    """create title blocks for run

    Args:
        bw_fname (str): name of brainweb file (MRI)
        MRE_fname (str): name of MRE file

    Returns:
        title, keyword file name, control file name
    """
    BW_short = bw_fname[0:9]
    MRE_short = "MRE" + MRE_fname[0:3]

    title = BW_short + "_" + MRE_short
    k_title = BW_short + "_" + MRE_short + ".k"
    control_title = "control_bw.k"

    return title, k_title, control_title


def local_fpaths():
    paths = {
        "BW_MRI_path": "T:/LSDYNA/brainwebworkflow/brainwebmri/",
        "BW_FE_path": "T:/LSDYNA/brainwebworkflow/brainwebFE/",
        "CTRL_FE_path": "T:/LSDYNA/brainwebworkflow/control_k/",
        "MRE_path": "T:/LSDYNA/brainwebworkflow/MRE134/",
        "img_path": "T:/LSDYNA/brainwebworkflow/coregistration_verify/",
        "icbm_path": "T:/LSDYNA/brainwebworkflow/ICBM/mni_icbm152_t1_tal_nlin_sym_09a.nii",
        "trans_path": "T:/LSDYNA/brainwebworkflow/transforms/",
    }

    return paths


def discovery_fpaths():
    paths = {
        "BW_MRI_path": "/work/muftu/Turner/brainwebworkflow/brainwebMRI/",
        "BW_FE_path": "/work/muftu/Turner/brainwebworkflow/brainwebFE/",
        "CTRL_FE_path": "/work/muftu/Turner/brainwebworkflow/control_k/",
        "MRE_path": "/work/muftu/Turner/brainwebworkflow/MRE134/",
        "img_path": "/work/muftu/Turner/brainwebworkflow/coregistration_verify/",
        "icbm_path": "/work/muftu/Turner/brainwebworkflow/ICBM/mni_icbm152_t1_tal_nlin_sym_09a.nii",
        "trans_path": "/work/muftu/Turner/brainwebworkflow/transforms/",
    }
    return paths


def run_file_names():
    bw_names = {
        "subject04": "BrainWeb_Subject04_updated.k",
        "subject05": "BrainWeb_Subject05_updated.k",
        "subject06": "BrainWeb_Subject06_updated.k",
        "subject18": "BrainWeb_Subject18_updated.k",
        "subject20": "BrainWeb_Subject20_updated.k",
        "subject38": "BrainWeb_Subject38_updated.k",
        "subject41": "BrainWeb_Subject41_updated.k",
        "subject42": "BrainWeb_Subject42_updated.k",
        "subject43": "BrainWeb_Subject43_updated.k",
        "subject44": "BrainWeb_Subject44_updated.k",
        "subject45": "BrainWeb_Subject45_updated.k",
        "subject46": "BrainWeb_Subject46_updated.k",
        "subject47": "BrainWeb_Subject47_updated.k",
        "subject48": "BrainWeb_Subject48_updated.k",
        "subject49": "BrainWeb_Subject49_updated.k",
        "subject50": "BrainWeb_Subject50_updated.k",
        "subject51": "BrainWeb_Subject51_updated.k",
        "subject52": "BrainWeb_Subject52_updated.k",
        "subject53": "BrainWeb_Subject53_updated.k",
        "subject54": "BrainWeb_Subject54_updated.k",
    }

    bw_mri_names = {
        "subject04": "subject04_crisp_v.mnc",
        "subject05": "subject05_crisp_v.mnc",
        "subject06": "subject06_crisp_v.mnc",
        "subject18": "subject18_crisp_v.mnc",
        "subject20": "subject20_crisp_v.mnc",
        "subject38": "subject38_crisp_v.mnc",
        "subject41": "subject41_crisp_v.mnc",
        "subject42": "subject42_crisp_v.mnc",
        "subject43": "subject43_crisp_v.mnc",
        "subject44": "subject44_crisp_v.mnc",
        "subject45": "subject45_crisp_v.mnc",
        "subject46": "subject46_crisp_v.mnc",
        "subject47": "subject47_crisp_v.mnc",
        "subject48": "subject48_crisp_v.mnc",
        "subject49": "subject49_crisp_v.mnc",
        "subject50": "subject50_crisp_v.mnc",
        "subject51": "subject51_crisp_v.mnc",
        "subject52": "subject52_crisp_v.mnc",
        "subject53": "subject53_crisp_v.mnc",
        "subject54": "subject54_crisp_v.mnc",
    }

    transforms = {
        "subject04": "0.nii",
        "subject05": "1.nii",
        "subject06": "2.nii",
        "subject18": "3.nii",
        "subject20": "4.nii",
        "subject38": "5.nii",
        "subject41": "6.nii",
        "subject42": "7.nii",
        "subject43": "8.nii",
        "subject44": "9.nii",
        "subject45": "10.nii",
        "subject46": "11.nii",
        "subject47": "12.nii",
        "subject48": "13.nii",
        "subject49": "14.nii",
        "subject50": "15.nii",
        "subject51": "16.nii",
        "subject52": "17.nii",
        "subject53": "18.nii",
        "subject54": "19.nii",
    }

    return bw_names, bw_mri_names, transforms


def COM_align(
    fixed: np.ndarray,
    moving: np.ndarray,
    fixed_mask: np.ndarray = None,
    moving_mask: np.ndarray = None,
):
    # Calculate Centers of Mass
    if fixed_mask is not None:
        COM_fixed = np.mean(fixed_mask, axis=0)
    else:
        COM_fixed = np.mean(fixed, axis=0)

    if moving_mask is not None:
        COM_moving = np.mean(moving_mask, axis=0)
    else:
        COM_moving = np.mean(moving, axis=0)

    # calculate translation and create transformation matrix
    offset = COM_fixed - COM_moving

    print(f"{datetime.now()}\t\tBrain COM offset by {offset}, correcting...")

    transform = np.array(
        [
            [0.8, 0, 0, offset[0]],
            [0, 0.8, 0, offset[1]],
            [0, 0, 0.8, offset[2]],
            [0, 0, 0, 1],
        ]
    )

    moving_augment = np.hstack((moving, np.ones((moving.shape[0], 1))))

    moving_transformed = moving_augment @ transform.T

    return moving_transformed[:, 0:3]


def find_keyword_includes(
    fpath: str, delimeter: str = "/"
) -> Union[str, list]:
    with open(fpath, "r") as f:
        ln = f.readlines()

    includes = []

    # find all includes, check for comments
    for idx, line in enumerate(ln):
        if "*INCLUDE" in line:
            if ln[idx + 1].startswith("$#"):
                includes.append(ln[idx + 2].split(delimeter)[-1])
            else:
                includes.append(ln[idx + 1].split(delimeter)[-1])

    # return appropriate values
    if len(includes) == 0:
        raise ValueError("No include statement found in file")
    elif len(includes) == 1:
        return includes[0]
    else:
        return includes


def clean_dir(path: str):
    # clean output directory

    spool_pattern = re.compile(r"spooles.res.\d{5}")
    mes_pattern = re.compile(r"mes\d{4}")

    num_cleaned = 0
    for filename in os.listdir(path):
        # Check if the filename matches the pattern
        if mes_pattern.match(filename):
            # Construct the full file path
            file_path = os.path.join(path, filename)
            # Delete the file
            os.remove(file_path)
            num_cleaned += 1
        elif spool_pattern.match(filename):
            # Construct the full file path
            file_path = os.path.join(path, filename)
            # Delete the file
            os.remove(file_path)
            num_cleaned += 1

    print(f"Cleaned {num_cleaned} files...")


def get_spacing(lims: dict, dims: tuple) -> tuple:
    maxes = lims["max"]

    mins = lims["min"]

    delta = maxes - mins

    spacing_x = delta[0] / dims[0]

    spacing_y = delta[1] / dims[1]

    spacing_z = delta[2] / dims[2]

    spc = (spacing_x, spacing_y, spacing_z)

    return spc


def parse_eigout(fpath: str, nfreqs: int = 25) -> np.ndarray:
    """Reads a LS-DYNA eigout file and returns the nfreqs lowest eingenvalue frequencies

    Args:
        fpath (str): Filepath to eigout file
        nfreqs (int, optional): Number of frequencies to return. Defaults to 25.

    Returns:
        np.ndarray: array of shape (nfreqs,) wit the lowest nfreqs eigenvalues
    """

    table_stem = (
        "        MODE    EIGENVALUE       RADIANS        CYCLES        PERIOD"
    )

    with open(fpath, "r") as f:
        lines = f.readlines()
        f.close()

    for idx, line in enumerate(lines):
        if table_stem in line:
            start_idx = idx + 2
            end_idx = start_idx + nfreqs
            break
    frequencies = []
    for n in range(start_idx, end_idx):
        linevals = lines[n].split()

        if len(linevals) >= 4:
            frequencies.append(float(linevals[3]))
        else:
            for i in range(n, end_idx):
                frequencies.append(0.0)
            break

    return np.array(frequencies)


def ants_affine(img: ants.core.ants_image.ANTsImage):
    # Extract image properties
    spacing = np.array(img.spacing)  # (sx, sy, sz)
    direction = np.array(img.direction).reshape((3, 3))  # 3x3 rotation matrix
    origin = np.array(img.origin)  # (ox, oy, oz)

    # Compute affine matrix
    affine_matrix = np.eye(4)  # Initialize as identity
    affine_matrix[:3, :3] = direction * spacing  # Scale direction by spacing
    affine_matrix[:3, 3] = origin  # Set translation

    return affine_matrix


if __name__ == "__main__":
    eigout_path = "kmeans/eigout"

    freqs = parse_eigout(eigout_path)
    print(freqs)
