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
import meshio
from ..models.femodel import FEModel
import numpy as np
from .calculate_prony import calculate_prony
from datetime import datetime
from typing import Union, List, Tuple
import os


def _entry_to_list(entry):
    """Check if an input is a single item and convert to a list if so.

    Args:
        entry (any): function input argument for parent function

    Returns:
        list: list form of input
    """

    if entry is not None and isinstance(entry, (str, ANTsImage)):
        return [entry]
    return entry


def _ensure_image(img: Union[str, ANTsImage]):
    """Ensure a correctly loaded Antsimage from an input which may
    contain an image or filepath string, or tuple of image or filepath strings

    Args:
        entry (str or AntsImage): Entry to check

    Raises:
        ValueError: File does not exist
        ValueError: entry is not str or AntsImage
    """
    if isinstance(img, str):
        if not os.path.exists(img):
            raise ValueError(f"No MRE image found at {img}")
        return image_read(img)
    elif isinstance(img, ANTsImage):
        return img
    else:
        raise ValueError("MRE_geom must include only str or AntsImage")


def _ensure_tuple(tup: Tuple[Union[str, ANTsImage]]):
    return tuple(_ensure_image(img) for img in tup)


def coregister_MRE_images(
    segmented_geom: Union[str, ANTsImage],
    target_label: int = 4,
    segmented_mask: Union[str, ANTsImage] = None,
    MRE_geom: List[Union[str, ANTsImage]] = None,
    MRE_mask: Union[str, ANTsImage] = None,
    MRE_to_transform: List[Tuple[Union[str, ANTsImage]]] = None,
    imgout: str = None,
    type_of_transform: str = "Affine",
):
    """Coregister MRE geometry image to segmented geometry image, and transform corresponding MRE images.

    Args:
        segmented_geom (Union[str, ANTsImage]): Segmented geometry image used for mesh creation
        target_label (int, optional): Label for ROI on geometry image. Used to create geometry mask if one is not provided. Defaults to 4.
        segmented_mask (Union[str,ANTsImage], optional): Binary mask for ROI on geometry image. Defaults to None.
        MRE_geom (List[Union[str, ANTsImage]], optional): List of geometry images associated with MRE at different frequencies. Defaults to None.
        MRE_mask (Union[str, ANTsImage], optional): Binary mask associated with ROI in MRE images. Defaults to None.
        MRE_to_transform (List[Tuple[Union[str, ANTsImage]]], optional): List of tuples, each tuple containing MRE images associated with one of the MRE_geom images provided. Defaults to None.
        imgout (str, optional): Filepath to save validation images of the transformations applied. Defaults to None.
        type_of_transform (str, optional): Type of transform, see Antspy registration documentation for guidance. Defaults to "Affine".

    Raises:
        ValueError: No MRE/segmented geometry image provided
        FileNotFoundError: MRE/segmented geometry files not found
        TypeError: Image files are not a filepath string or AntsImage
        ValueError: Transform could not be resolved
        ValueError: MRE_geom and MRE_to_Transform are different lengths
        ValueError: imgout path could not be found

    Returns:
        List: list of transformations
        List[tuple]: list of transformed image tuples
    """
    if segmented_geom is None:
        raise ValueError("Segmented geometry image is required")

    if MRE_geom is None:
        raise ValueError("MRE geometry image is required")

    MRE_geom = _entry_to_list(MRE_geom)
    MRE_to_transform = _entry_to_list(MRE_to_transform)

    # Load segmented geometry
    if isinstance(segmented_geom, str):
        if not os.path.exists(segmented_geom):
            raise FileNotFoundError(
                f"Geometry mask file not found: {segmented_geom}"
            )
        segmented_geom = image_read(segmented_geom)
    elif not isinstance(segmented_geom, ANTsImage):
        raise TypeError(
            "geom_mask must be either a filepath string or ANTsImage object"
        )

    # Load segmented mask or threshold image
    if segmented_mask is None:
        segmented_mask = threshold_image(
            segmented_geom,
            low_thresh=target_label - 0.01,
            high_thresh=target_label + 1.01,
        )
    elif isinstance(segmented_mask, str):
        if not os.path.exists(segmented_mask):
            raise FileNotFoundError(
                f"Geometry image file not found: {segmented_mask}"
            )
        segmented_geom = image_read(segmented_geom)
    elif not isinstance(segmented_geom, ANTsImage):
        raise TypeError(
            "geom must be either a filepath string or ANTsImage object"
        )

    # load MRE geometries
    MRE_geom_imgs = []
    for img in MRE_geom:
        MRE_geom_imgs.append(_ensure_image(img))

    # load MRE mask
    if MRE_mask is not None and isinstance(MRE_mask, str):
        if not os.path.exists(MRE_mask):
            raise FileNotFoundError(
                f"Geometry image file not found: {MRE_mask}"
            )
        MRE_mask = image_read(MRE_mask)
    elif not isinstance(segmented_geom, ANTsImage):
        raise TypeError(
            "geom must be either a filepath string or ANTsImage object"
        )

    # load images to transform
    MRE_transform_imgs = []
    for imgs in MRE_to_transform:
        MRE_transform_imgs.append(_ensure_tuple(imgs))

    if len(MRE_transform_imgs) != len(MRE_geom_imgs):
        raise ValueError(
            f"{len(MRE_geom_imgs)} geometry images were provided but {len(MRE_transform_imgs)} transform tuples were provided"
        )

    # create each transform and image
    transformations = []
    transformed_images = []
    for idx, (geom, img_tuple) in enumerate(zip(MRE_geom, MRE_transform_imgs)):
        # resample geometry to MRE image
        geom_resample = resample_image_to_target(segmented_geom, geom)
        mask_resample = resample_image_to_target(segmented_mask, geom)

        # calculate transform
        try:
            if MRE_mask is not None:
                tx = registration(
                    fixed=geom_resample,
                    moving=geom,
                    type_of_transform=type_of_transform,
                    mask=mask_resample,
                    moving_mask=MRE_mask,
                )
            else:
                tx = registration(
                    fixed=geom_resample,
                    moving=geom,
                    type_of_transform=type_of_transform,
                    mask=mask_resample,
                )
        except ValueError:
            print(f"transformation failed on image {idx}")

        transformations.append(tx["fwdtransforms"])

        # apply transforms
        transformed_tuple = tuple(
            apply_transforms(
                fixed=geom_resample,
                moving=img,
                transformlist=tx["fwdtransforms"],
            )
            for img in img_tuple
        )
        transformed_images.append(transformed_tuple)

        if imgout is not None:
            if not os.path.exists(imgout):
                raise ValueError("imgout directory does not exist")
            else:
                base = f"{imgout + 'MRE{idx}_coreg.jpg'}"
                plot(
                    segmented_geom,
                    overlay=tx["warpedmovout"],
                    overlay_cmap="Dark2",
                    overlay_alpha=0.8,
                    filename=base,
                    axis=0,
                )

    # return single dict or list of dicts
    if len(transformed_images) == 0:
        raise ValueError("No results generated from MRE coregistration")
    elif len(transformed_images) == 1:
        return transformations[0], transformed_images[0]
    else:
        return transformations, transformed_images


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
        print(
            f"Segment {i}: SS={SS_mean:.3f}, DR={DR_mean:.3f}, Ginf={Ginf:.3f}, G1={G1:.3f}, tau={tau:.3f}"
        )

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
    coreg_results = coregister_MRE_images(
        segmented_geom=geom, gp_list=SS_list, gpp_list=DR_list, imgout=imgout
    )
    all_results = []
    for result in coreg_results:
        ROI_means = segment_MRE_regions(result["gp"], result["gpp"], n_segs)
        all_results.append(ROI_means)
    return coreg_results, all_results


def meshio_to_femodel(
    mesh: meshio.Mesh,
    title: str = "",
    source: str = "",
    default_part_id: int = 1,
) -> FEModel:
    """
    Convert a meshio.Mesh object into a custom FEModel object.

    Args:
        mesh: meshio.Mesh object
        title: Metadata title for FEModel
        source: Metadata source for FEModel
        default_part_id: Default part ID for all elements

    Returns:
        FEModel instance with custom nodes and elements
    """
    femodel = FEModel(title=title, source=source)

    # Add nodes
    for node_id, (x, y, z) in enumerate(mesh.points, start=1):
        femodel.add_node(node_id, x, y, z)

    # Handle only one type of element for now (ex. "tetra")
    supported_keys = ["tetra", "triangle", "hexahedron", "quad"]
    found = False
    for key in supported_keys:
        if key in mesh.cells_dict:
            elements = mesh.cells_dict[key]
            for elem_id, node_ids in enumerate(elements, start=1):
                # FIXED: + 1 offset since meshio is zero indexed but FEModel is one indexed
                femodel.add_element(
                    elem_id, [i + 1 for i in node_ids], default_part_id
                )
            found = True
            break

    if not found:
        raise ValueError(
            f"No supported cell types found in mesh. Supported: {supported_keys}"
        )

    return femodel
