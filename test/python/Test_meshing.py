import pytest

import numpy as np

import os

import meshio

from MRI2FE import mesh_from_nifti, nifti_to_inr, FEModel, model_from_meshio


class TestMeshing:
    def test_file_conversion(self):
        root_dir = os.getcwd()

        path = os.path.join(
            root_dir, "test", "test_data", "test_concentric_spheres.nii"
        )

        outpath = nifti_to_inr(path)

        assert os.path.exists(outpath)

    def test_file_raises(self):
        root_dir = os.getcwd()

        path = os.path.join(root_dir, "test", "test_data", "not_real_file.nii")
        with pytest.raises(ValueError):
            outpath = nifti_to_inr(path)

    def test_mesh_creation(self):
        root_dir = os.getcwd()

        path = os.path.join(
            root_dir, "test", "test_data", "test_concentric_spheres.nii"
        )

        mesh = mesh_from_nifti(path, optimize=True, facetSize=3, cellSize=2)

        assert mesh is not None

        assert mesh.points.shape[0] > 0
        assert mesh.points.shape[1] > 0

        out_path = os.path.join(
            root_dir, "test", "test_data", "test_mesh_out.mesh"
        )

        mesh.write(out_path)

        centered_points = mesh.points - np.array([32.0, 32.0, 32.0])

        dist = np.linalg.norm(centered_points, axis=1)

        ref_dist = 33 * np.ones_like(dist)

        np.testing.assert_array_less(dist, ref_dist)

    def test_femodel_from_meshio(self):
        root_dir = os.getcwd()

        inpath = os.path.join(
            root_dir, "test", "test_data", "test_mesh_out.mesh"
        )

        mesh = meshio.read(inpath)
        print(mesh.points)
        print(mesh.cells_dict["tetra"])

        mdl: FEModel = model_from_meshio(mesh, title="test", source="test")

        print(mdl.get_node_table().shape)
        print(mdl.get_element_table().shape)

        assert mdl.get_node_table().shape[0] > 1
