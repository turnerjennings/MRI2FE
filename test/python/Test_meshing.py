import pytest

import numpy as np

import os

import meshio

from ants import image_read

from MRI2FE import mesh_from_nifti, nifti_to_inr, FEModel, spatial_map


class TestMeshing:
    def test_file_conversion(self):
        root_dir = os.getcwd()

        path = os.path.join(
            root_dir, "test", "test_data", "test_concentric_spheres.nii"
        )

        img = image_read(path)

        outpath = nifti_to_inr(img)

        assert os.path.exists(outpath)

    def test_file_raises(self):
        root_dir = os.getcwd()

        path = os.path.join(root_dir, "test", "test_data", "not_real_file.nii")
        with pytest.raises(ValueError):
            outpath = mesh_from_nifti(path)

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

        # create spatial map from original image
        img_map = spatial_map(image_read(path))

        # img_map = img_map[img_map[:,3] > 0,:]

        img_map_min_max = [
            np.min(img_map[:, :3], axis=0),
            np.max(img_map[:, :3], axis=0),
        ]

        mesh_min_max = [
            np.min(mesh.points, axis=0),
            np.max(mesh.points, axis=0),
        ]

        print(image_read(path))
        print("nifti min max:")
        print(img_map_min_max)
        print("mesh min max:")
        print(mesh_min_max)

        np.testing.assert_array_less(img_map_min_max[0], mesh_min_max[0])
        np.testing.assert_array_less(mesh_min_max[1], img_map_min_max[1])

    def test_femodel_from_meshio(self):
        root_dir = os.getcwd()

        inpath = os.path.join(
            root_dir, "test", "test_data", "test_mesh_out.mesh"
        )

        mesh: meshio.Mesh = meshio.read(inpath)

        for idx, item in enumerate(mesh.cells):
            if item.type == "tetra":
                shp = item.data.shape

        mdl = FEModel(title="test", source="test")

        mdl.from_meshio(mesh)


        assert mdl.node_table.shape[0] == mesh.points.shape[0]

        assert mdl.element_table.shape[0] == shp[0]

        # check for zero nodes and node range
        assert np.min(mdl.node_table[:, 0]) > 0

        assert np.min(mdl.element_table[:, 2:]) == np.min(mdl.node_table[:, 0])

        assert np.max(mdl.element_table[:, 2:]) == np.max(mdl.node_table[:, 0])

        # check mesh regions
        shape = (64, 64, 64)
        radii = [5, 10, 15, 20]  # Radii for spheres
        center = np.array(shape) // 2
        origin = (4.0, -1.0, 2.0)

        mdl.update_centroids()

        center_dist_mesh = np.linalg.norm(
            mdl.centroid_table - (center + origin), axis=1
        )

        for i in range(1, 5):
            print(i)
            points_in_region = mdl.element_table[:, 1] == i
            print(mdl.centroid_table[points_in_region, :])
            print(center + origin)
            print(mdl.centroid_table[points_in_region, :] - (center + origin))
            print(center_dist_mesh[points_in_region])

            np.testing.assert_array_less(
                center_dist_mesh[points_in_region], radii[i - 1] + 0.1
            )
