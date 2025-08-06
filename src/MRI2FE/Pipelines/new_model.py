from typing import Any, List, Literal, Optional, Tuple

from ants import image_read

from ..generate_mesh import mesh_from_nifti
from ..models.femodel import FEModel
from ..MRE.calculate_prony import calculate_prony
from ..MRE.MRE_coregistration import coregister_MRE_images, segment_MRE_regions
from ..MRE.MRE_mapping import map_MRE_to_mesh


class FEModelbuilder:
    def __init__(self, title: str = "", source: str = ""):
        """Initialize the FEModel object to store model data.

        Args:
            title (str, optional): Optional title for the model which will be written to output solver decks. Defaults to "".
            source (str, optional): Optional source folder for model for internal tracking. Defaults to "".
        """
        self.model = FEModel(title=title, source=source)

    def mesh(
        self,
        img_path: str,
        img_labels: Optional[List[str]] = None,
        optimize: bool = False,
        **kwargs,
    ):
        """Generate a tetrahedral mesh from labeled MRI data and store in the FEModel object

        Args:
            img_path (str): Path to segmented, labeled MRI image.
            img_labels (List[str], optional): Optional labels providing names for each region of the labeled image. Defaults to None.
            optimize (bool, optional): Whether to perform post-process optimization on the tetrahedral mesh.  Increases quality and run time. Defaults to False.
        """
        self.labeled_geom = image_read(img_path)
        self.geom_labels = img_labels

        msh = mesh_from_nifti(filepath=img_path, optimize=optimize, **kwargs)

        self.model.from_meshio(msh, region_names=img_labels)

        return self

    def map_mre(
        self,
        target_label: int = 4,
        MRE_type: Literal[
            "stiffness_damping", "complex_shear"
        ] = "stiffness_damping",
        MRE_geom: Optional[List[str | Any]] = None,
        MRE_mask: Optional[str | Any] = None,
        MRE_frequency: Optional[List[float]] = None,
        MRE_to_transform: Optional[List[Tuple[str | Any]]] = None,
        n_segs: int = 5,
        **kwargs,
    ):
        """Calculate material model coefficients and map MRE material assignments onto an ROI on the mesh.

        Args:
            target_label (int, optional): Target integer label on the labeled image to map MRE material properties to. Defaults to 4.
            MRE_type ("stiffness_damping" or "complex_shear", optional): Specify whether MRE files are provided as shear stiffness and damping ratio or as storage and loss moduli. Defaults to "stiffness_damping".
            MRE_geom (List[str  |  Any], optional): List of images or paths to images for MRE geometries at each frequency. Defaults to None.
            MRE_mask (str | Any, optional): ROI mask for MRE geometry. Defaults to None.
            MRE_frequency (List[float], optional): List of frequencies for each MRE geometry image. Defaults to None.
            MRE_to_transform (List[Tuple[str  |  Any]], optional): List of tuples of strings or MRE images.  Each tuple represents either shear/damping or storage/loss moduli at a frequency given in MRE_frequency. Defaults to None.
            n_segs (int, optional): Number of segments to discretize the MRE material properties into. Defaults to 5.

        """
        _, transformed = coregister_MRE_images(
            segmented_geom=self.labeled_geom,
            target_label=target_label,
            MRE_geom=MRE_geom,
            MRE_mask=MRE_mask,
            MRE_to_transform=MRE_to_transform,
            **kwargs,
        )
        # handle edge case if only one MRE frequency is used
        if isinstance(transformed, tuple):
            transformed = [transformed]

        self.transformed_mre = transformed

        labels, region_avgs = segment_MRE_regions(
            img_list=transformed, n_segs=n_segs
        )

        regions_props = []
        for i in range(len(region_avgs)):
            if MRE_type == "stiffness_damping":
                regions_props.append(
                    calculate_prony(
                        mu=region_avgs["1"][i],
                        xi=region_avgs["2"][i],
                        w=MRE_frequency,
                    )
                )

            elif MRE_type == "complex_shear":
                regions_props.append(
                    calculate_prony(
                        gp=region_avgs["1"][i],
                        gpp=region_avgs["2"][i],
                        w=MRE_frequency,
                    )
                )

        if self.geom_labels is not None:
            self.model = map_MRE_to_mesh(
                self.model,
                label_img=labels,
                region_properties=regions_props,
                target_region_id=target_label,
                region_prefix=self.geom_labels[target_label - 1],
            )
        else:
            self.model = map_MRE_to_mesh(
                self.model,
                label_img=labels,
                region_properties=regions_props,
                target_region_id=target_label,
            )

        return self

    def write(self, fpath: str, type: Literal["lsdyna"] = "lsdyna"):
        """Write model to output solver deck

        Args:
            fpath (str): File path to save output to.
            type ("lsdyna", optional): Output type to be saved.  Currently, only LS-DYNA is supported. Defaults to "lsdyna".

        """
        if type == "lsdyna":
            self.model.write_lsdyna(fpath)

        return self

    def build(self):
        """Return generated FEModel"""
        return self.model
