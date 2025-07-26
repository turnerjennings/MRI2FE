from MRI2FE.Pipelines.new_model import FEModelbuilder
import os
from ants import image_read
import pytest


class Test_FEModelbuilder:
    @pytest.fixture(autouse=True)
    def define_fpaths(self):
        self.labeled_geom_path = os.path.join(
            "test", "test_data", "test_concentric_spheres.nii"
        )

        self.MRE_geom_path = os.path.join(
            "test", "test_data", "test_concentric_spheres.nii"
        )

        self.MRE_ss_path = os.path.join(
            "test", "test_data", "test_stiffness.nii"
        )

        self.MRE_dr_path = os.path.join(
            "test", "test_data", "test_damping_ratio.nii"
        )

    def test_pipeline_mesh(self):
        test_mdl = (
            FEModelbuilder()
            .mesh(
                img_path=self.labeled_geom_path,
                img_labels=["one", "two", "three", "four"],
            )
            .build()
        )

        assert test_mdl.metadata["num_nodes"] > 0
        assert test_mdl.metadata["num_elements"] > 0
        assert test_mdl.node_table.size > 0
        assert test_mdl.element_table.size > 0

        assert len(test_mdl.part_info) > 0

    def test_pipeline_mre(self):
        freq = [30]

        MRE_to_transform = [(self.MRE_ss_path, self.MRE_dr_path)]

        test_mdl = (
            FEModelbuilder()
            .mesh(
                img_path=self.labeled_geom_path,
                img_labels=["one", "two", "three", "four"],
            )
            .map_mre(
                target_label=1,
                MRE_type="stiffness_damping",
                MRE_geom=self.MRE_geom_path,
                MRE_mask=None,
                MRE_frequency=freq,
                MRE_to_transform=MRE_to_transform,
            )
            .build()
        )

        print(test_mdl.material_info)
        print(test_mdl.part_info)

        assert len(test_mdl.material_info) == 1
        assert len(test_mdl.part_info) == 4 + 1

    def test_pipeline_write(self, tmp_path):
        mu = image_read(self.MRE_ss_path)
        xi = image_read(self.MRE_dr_path)

        freq = [30]

        MRE_to_transform = [(mu, xi)]

        _ = (
            FEModelbuilder()
            .mesh(
                img_path=self.labeled_geom_path,
                img_labels=["one", "two", "three", "four"],
            )
            .map_mre(
                target_label=1,
                MRE_type="stiffness_damping",
                MRE_geom=self.MRE_geom_path,
                MRE_mask=None,
                MRE_frequency=freq,
                MRE_to_transform=MRE_to_transform,
            )
            .write(os.path.join(tmp_path, "test_out.k"))
            .build()
        )

        assert os.path.exists(os.path.join(tmp_path, "test_out.k"))
