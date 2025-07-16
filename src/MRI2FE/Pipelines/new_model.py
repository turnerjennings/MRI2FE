from typing import Literal, Dict, Any
from ants.core.ants_image import ANTsImage
from ants import image_read
from ..MRE.MRE_coregistration import coregister_MRE_images, segment_MRE_regions
from ..MRE.calculate_prony import calculate_prony
from ..output.k_file_operations import parse_k_file
from ..MRE.MRE_mapping import map_MRE_to_mesh
from ..FEModel.femodel import FEModel
from ..utilities import element_centroids
import numpy as np


class NewModel:
    def mre_mri_registration(
        self,
        mri_geometry_path: str,
        mre_type: Literal["complex_shear", "stiffness_damping"],
        coregister_mre: bool = True,
        n_segments: int = 5,
        complex_modulus_path: str = None,
        stiffness_path: str = None,
        damping_ratio_path: str = None,
        output_viz_path: str = None,
    ):
        """Register MRE images to MRI geometry and optionally segment them.

        Args:
            mri_geometry_path: Path to MRI geometry file (.nii)
            mre_type: Format of MRE data ("complex_shear" or "stiffness_damping")
            coregister_mre: Whether to perform coregistration
            n_segments: Number of segments for MRE segmentation
            complex_modulus_path: Path to complex shear modulus file (if mre_type is "complex_shear")
            stiffness_path: Path to stiffness file (if mre_type is "stiffness_damping")
            damping_ratio_path: Path to damping ratio file (if mre_type is "stiffness_damping")
            output_viz_path: Optional path to save visualization images

        Returns:
            dict: Registered and/or segmented images and transformation data
        """

        geom_image = image_read(mri_geometry_path)

        if not coregister_mre:
            return {"geometry": geom_image}

        if mre_type == "complex_shear":
            if not complex_modulus_path:
                raise ValueError(
                    "complex_modulus_path required for complex_shear MRE type"
                )

            result = coregister_MRE_images(
                geom=geom_image,
                geom_mask=None,
                gp_list=complex_modulus_path,
                gpp_list=complex_modulus_path,
                mu_list=None,
                xi_list=None,
                imgout=output_viz_path,
            )

        elif mre_type == "stiffness_damping":
            if not (stiffness_path and damping_ratio_path):
                raise ValueError(
                    "Both stiffness_path and damping_ratio_path required for stiffness_damping MRE type"
                )

            result = coregister_MRE_images(
                geom=geom_image,
                geom_mask=None,
                gp_list=None,
                gpp_list=None,
                mu_list=stiffness_path,
                xi_list=damping_ratio_path,
                imgout=output_viz_path,
            )

        return result

    def segment_MRE_images(
        self,
        registered_images: Dict[str, Any],
        mre_type: Literal["complex_shear", "stiffness_damping"],
        n_segments: int = 5,
    ) -> Dict[str, Any]:
        """Segment the registered MRE images into regions.

        Args:
            registered_images: Dictionary containing registered images from mre_mri_registration
            mre_type: Format of MRE data ("complex_shear" or "stiffness_damping")
            n_segments: Number of segments to create (default: 5)

        Returns:
            dict: Contains segmentation results and Prony parameters for each region
        """
        if not registered_images:
            raise ValueError("registered_images is required")

        if mre_type == "complex_shear":
            if "gp" not in registered_images or "gpp" not in registered_images:
                raise ValueError(
                    "registered_images must contain 'gp' and 'gpp' for complex_shear type"
                )

            return segment_MRE_regions(
                SS_img=registered_images["gp"],
                DR_img=registered_images["gpp"],
                n_segs=n_segments,
            )

        elif mre_type == "stiffness_damping":
            if "mu" not in registered_images or "xi" not in registered_images:
                raise ValueError(
                    "registered_images must contain 'mu' and 'xi' for stiffness_damping type"
                )

            return segment_MRE_regions(
                SS_img=registered_images["mu"],
                DR_img=registered_images["xi"],
                n_segs=n_segments,
            )

    def calculate_material_constants(
        self,
        registered_images: Dict[str, ANTsImage],
        mre_type: Literal["complex_shear", "stiffness_damping"],
        mre_frequency: float = 50.0,
    ) -> Dict[str, list]:
        """Calculate material constants (Prony series parameters) using the Prony series model.

        Args:
            registered_images: Dictionary containing registered images from mre_mri_registration
            mre_type: Format of MRE data ("complex_shear" or "stiffness_damping")
            mre_frequency: MRE vibration frequency in Hz (default: 50.0)

        Returns:
            dict: Contains material constants:
                - Ginf: Long-term shear modulus
                - G1: Short-term shear modulus
                - Tau: Relaxation time constant
        """
        if not registered_images:
            raise ValueError("registered_images is required")

        if mre_type == "complex_shear":
            if "gp" not in registered_images or "gpp" not in registered_images:
                raise ValueError(
                    "registered_images must contain 'gp' and 'gpp' for complex_shear type"
                )

            gp_array = registered_images["gp"].numpy()
            gpp_array = registered_images["gpp"].numpy()

            nonzero_mask = (gp_array > 0) & (gpp_array > 0)
            gp_values = gp_array[nonzero_mask]
            gpp_values = gpp_array[nonzero_mask]

            if len(gp_values) == 0 or len(gpp_values) == 0:
                raise ValueError(
                    "No valid non-zero values found in the images"
                )

            w_array = np.full_like(gp_values, mre_frequency)

            Ginf, G1, tau = calculate_prony(
                gp=gp_values, gpp=gpp_values, w=w_array, tol=1.0
            )

        elif mre_type == "stiffness_damping":
            if "mu" not in registered_images or "xi" not in registered_images:
                raise ValueError(
                    "registered_images must contain 'mu' and 'xi' for stiffness_damping type"
                )

            mu_array = registered_images["mu"].numpy()
            xi_array = registered_images["xi"].numpy()

            nonzero_mask = (mu_array > 0) & (xi_array > 0)
            mu_values = mu_array[nonzero_mask]
            xi_values = xi_array[nonzero_mask]

            if len(mu_values) == 0 or len(xi_values) == 0:
                raise ValueError(
                    "No valid non-zero values found in the images"
                )

            mu_values = mu_values / 1000.0
            w_array = np.full_like(mu_values, mre_frequency)

            try:
                Ginf, G1, tau = calculate_prony(
                    mu=mu_values, xi=xi_values, w=w_array, tol=1.0
                )
            except ValueError:
                mu_mean = float(np.mean(mu_values))
                xi_mean = float(np.mean(xi_values))

                Ginf, G1, tau = calculate_prony(
                    mu=np.array([mu_mean]),
                    xi=np.array([xi_mean]),
                    w=np.array([mre_frequency]),
                    tol=1.0,
                )
            except Exception as e:
                print(f"Error calculating Prony series: {str(e)}")
                raise e

        return {"Ginf": Ginf, "G1": G1, "Tau": tau}

    def calculate_element_centroids(
        self, mesh_path: str
    ) -> Dict[str, np.ndarray]:
        """Calculate centroids for all elements in the FE mesh.

        Args:
            mesh_path: Path to the FE mesh file (.k format)

        Returns:
            dict: Contains element data:
                - element_table: Array of element connectivity
                - node_table: Array of node coordinates
                - centroids: Array of element centroids
        """
        if not mesh_path:
            raise ValueError("mesh_path is required")

        element_table, node_table = parse_k_file(mesh_path)

        n_elements = len(element_table)
        centroids = np.zeros((n_elements, 3))

        for i, element in enumerate(element_table):
            try:
                centroid = element_centroids(element, node_table)
                centroids[i] = centroid
            except Exception as e:
                raise e

        return {
            "element_table": element_table,
            "node_table": node_table,
            "centroids": centroids,
        }

    def assign_material_properties(
        self,
        material_constants: Dict[str, Any],
        mesh_data: Dict[str, np.ndarray],
        mre_map: ANTsImage,
        offset: int = 3,
        output_csv: str = None,
    ) -> FEModel:
        """Assign material properties to the FE mesh based on MRE data.

        Args:
            material_constants: Dictionary containing Prony series parameters from calculate_material_constants
            mesh_data: Dictionary containing mesh data from calculate_element_centroids
            mre_map: Segmented MRE image map
            offset: Optional filter for part IDs (default: 3)
            output_csv: Optional path to save CSV output files

        Returns:
            FEModel: Updated finite element model with assigned material properties
        """
        if not material_constants or not mesh_data:
            raise ValueError("material_constants and mesh_data are required")

        if not isinstance(mre_map, ANTsImage):
            raise ValueError("mre_map must be an ANTsImage")

        fe_model = FEModel()

        for node_data in mesh_data["node_table"]:
            node_id = int(node_data[0])
            x, y, z = node_data[1:]
            fe_model.add_node(node_id=node_id, x=x, y=y, z=z)

        updated_model = map_MRE_to_mesh(
            fe_model=fe_model,
            map=mre_map,
            elcentroids=mesh_data["centroids"],
            ect=mesh_data["element_table"],
            offset=offset,
            csvpath=output_csv,
        )

        return updated_model
