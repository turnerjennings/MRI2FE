from typing import Literal, Any
from ants import image_read
from ..MRE.MRE_coregistration import coregister_MRE_images, segment_MRE_regions
from ..MRE.calculate_prony import calculate_prony
from ..MRE.MRE_mapping import map_MRE_to_mesh
from ..models.femodel import FEModel
from ..generate_mesh import mesh_from_nifti
from typing import List, Tuple


class FEModelbuilder:
    def __init__(self):
        self.model = FEModel()

    def mesh(
        self,
        img_path: str,
        img_labels: List[str] = None,
        optimize: bool = False,
        **kwargs,
    ):
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
        MRE_geom: List[str | Any] = None,
        MRE_mask: str | Any = None,
        MRE_frequency: List[float] = None,
        MRE_to_transform: List[Tuple[str | Any]] = None,
        n_segs: int = 5,
        **kwargs,
    ):
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
        if type == "lsdyna":
            self.model.write_lsdyna(fpath)

        return self

    def build(self):
        return self.model
