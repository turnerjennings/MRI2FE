import pytest
import numpy as np

from ants import image_read, from_numpy

from MRI2FE.MRE import (
    coregister_MRE_images,
    segment_MRE_regions,
    map_MRE_to_mesh,
)
from MRI2FE import FEModel, spatial_map

import os

import meshio


class TestMRECoreg:
    @pytest.fixture(autouse=True)
    def problem_setup(self):
        self.geom_path = os.path.join(
            "test", "test_data", "test_concentric_spheres.nii"
        )

        self.MRE_geom_path = os.path.join(
            "test", "test_data", "test_concentric_spheres.nii"
        )

        self.gp_path = os.path.join("test", "test_data", "test_stiffness.nii")

        self.gpp_path = os.path.join(
            "test", "test_data", "test_damping_ratio.nii"
        )

        self.geom = image_read(self.geom_path)

        geom_data = self.geom.numpy()

        geom_data[geom_data > 0] = 1

        self.geom_mask = from_numpy(
            geom_data,
            origin=self.geom.origin,
            spacing=self.geom.spacing,
            direction=self.geom.direction,
        )

        self.MRE_geom = image_read(self.MRE_geom_path)

        geom_data = self.MRE_geom.numpy()

        geom_data[geom_data > 0] = 1

        self.MRE_mask = from_numpy(
            geom_data,
            origin=self.MRE_geom.origin,
            spacing=self.MRE_geom.spacing,
            direction=self.MRE_geom.direction,
        )

        self.gp = image_read(self.gp_path)
        self.gpp = image_read(self.gpp_path)

    def test_coreg_single(self):
        MRE_images = [(self.gp, self.gpp)]

        transforms, images = coregister_MRE_images(
            segmented_geom=self.geom,
            segmented_mask=self.geom_mask,
            MRE_geom=self.MRE_geom,
            MRE_mask=self.MRE_mask,
            MRE_to_transform=MRE_images,
        )

        assert len(transforms) == 1
        assert len(images) == 2

    def test_coreg_single_fromstring(self):
        MRE_images = [(self.gp_path, self.gpp_path)]

        transforms, images = coregister_MRE_images(
            segmented_geom=self.geom,
            segmented_mask=self.geom_mask,
            MRE_geom=self.MRE_geom,
            MRE_mask=self.MRE_mask,
            MRE_to_transform=MRE_images,
        )

        assert len(transforms) == 1
        assert len(images) == 2

    def test_coreg_multiple(self):
        MRE_images = [
            (self.gp, self.gpp),
            (self.gp, self.gpp),
            (self.gp, self.gpp),
            (self.gp, self.gpp),
        ]

        MRE_geoms = [
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
        ]

        transforms, images = coregister_MRE_images(
            segmented_geom=self.geom,
            segmented_mask=self.geom_mask,
            MRE_geom=MRE_geoms,
            MRE_mask=self.MRE_mask,
            MRE_to_transform=MRE_images,
        )

        assert len(transforms) == 4
        assert len(images) == 4

    def test_coreg_multiple_frompath(self):
        MRE_images = [
            (self.gp_path, self.gpp_path),
            (self.gp_path, self.gpp_path),
            (self.gp_path, self.gpp_path),
            (self.gp_path, self.gpp_path),
        ]

        MRE_geoms = [
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
        ]

        transforms, images = coregister_MRE_images(
            segmented_geom=self.geom,
            segmented_mask=self.geom_mask,
            MRE_geom=MRE_geoms,
            MRE_mask=self.MRE_mask,
            MRE_to_transform=MRE_images,
        )

        assert len(transforms) == 4
        assert len(images) == 4

    def test_raises(self):
        MRE_images = [
            (self.gp, self.gpp),
            (self.gp, self.gpp),
            (self.gp, self.gpp),
            (self.gp, self.gpp),
        ]

        MRE_geoms = [
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
            self.MRE_geom,
        ]

        # no geometry provided
        with pytest.raises(TypeError):
            transforms, images = coregister_MRE_images(
                segmented_mask=self.geom_mask,
                MRE_geom=MRE_geoms,
                MRE_mask=self.MRE_mask,
                MRE_to_transform=MRE_images,
            )

        # no MRE geometry provided
        with pytest.raises(ValueError):
            transforms, images = coregister_MRE_images(
                segmented_geom=self.geom,
                segmented_mask=self.geom_mask,
                MRE_mask=self.MRE_mask,
                MRE_to_transform=MRE_images,
            )

        # files not found
        with pytest.raises(FileNotFoundError):
            transforms, images = coregister_MRE_images(
                segmented_geom="nonexistent/file/path",
                segmented_mask=self.geom_mask,
                MRE_geom=MRE_geoms,
                MRE_mask=self.MRE_mask,
                MRE_to_transform=MRE_images,
            )

        MRE_geoms = [self.MRE_geom, self.MRE_geom, self.MRE_geom]

        # MRE and MRE_to_Transform different lengths

        with pytest.raises(ValueError):
            transforms, images = coregister_MRE_images(
                segmented_geom=self.geom,
                segmented_mask=self.geom_mask,
                MRE_geom=MRE_geoms,
                MRE_mask=self.MRE_mask,
                MRE_to_transform=MRE_images,
            )


class TestMRESegmentation:
    def test_seg_single(self):
        ss = image_read("test/test_data/test_stiffness.nii")

        dr = image_read("test/test_data/test_damping_ratio.nii")

        img_list = [(ss, dr)]

        test_img, test_dict = segment_MRE_regions(img_list=img_list, n_segs=5)

        print(test_dict["1"])
        print(test_dict["2"])

        assert len(test_dict["1"]) == 6
        assert len(test_dict["2"]) == 6

        for list in test_dict["1"]:
            assert len(list) == 1

        for list in test_dict["2"]:
            assert len(list) == 1

        assert test_img.max() == 5

    def test_seg_multiple(self):
        ss = image_read("test/test_data/test_stiffness.nii")

        dr = image_read("test/test_data/test_damping_ratio.nii")

        img_list = [(ss, dr), (ss, dr), (ss, dr)]

        test_img, test_dict = segment_MRE_regions(img_list=img_list, n_segs=5)

        assert len(test_dict["1"]) == 6
        assert len(test_dict["2"]) == 6

        for list in test_dict["1"]:
            assert len(list) == 3

        for list in test_dict["2"]:
            assert len(list) == 3

        assert test_img.max() == 5


class TestMREMapping:
    @pytest.fixture(autouse=True)
    def create_MRE_labels(self):
        # Parameters
        shape = (64, 64, 64)
        center = np.array(shape) // 2
        radii = [2, 3, 4]  # Radii for spheres
        labels = list(range(1, len(radii) + 1))  # Labels: 1, 2, 3, ...

        spacing = (1.0, 1.0, 1.0)
        origin = (4.0, -1.0, 2.0)

        # Create empty volume
        volume = np.zeros(shape, dtype=np.uint8)

        # Grid of coordinates
        zz, yy, xx = np.meshgrid(
            np.arange(shape[0]),
            np.arange(shape[1]),
            np.arange(shape[2]),
            indexing="ij",
        )
        dist = np.sqrt(
            (xx - center[2]) ** 2
            + (yy - center[1]) ** 2
            + (zz - center[0]) ** 2
        )

        # Assign labels to voxels inside each sphere
        for i, r in enumerate(radii):
            if i == 0:
                mask = dist <= r
            else:
                mask = (dist <= r) & (dist > radii[i - 1])
            volume[mask] = labels[i]

        # Save as NIfTI
        self.mre_labels = from_numpy(
            data=volume, origin=origin, spacing=spacing
        )
        self.obj_center = center
        self.obj_origin = origin
        self.mapping_radii = radii
        self.n_mre_regions = len(radii)

    @pytest.fixture(autouse=True)
    def load_mesh(self):
        mesh_path = os.path.join("test", "test_data", "test_mesh_out.mesh")

        mesh = meshio.read(mesh_path)

        self.mdl = FEModel(title="test_mesh", source=mesh_path)

        self.mdl.from_meshio(mesh)

        self.mdl.update_centroids()

    def test_check_mesh_inputs(self):
        min_max_mesh = [
            np.min(self.mdl.node_table[:, 1:], axis=0),
            np.max(self.mdl.node_table[:, 1:], axis=0),
        ]

        print(f"Mesh min/maxes: {min_max_mesh}")

        min_max_img = [
            self.mre_labels.origin,
            np.array(self.mre_labels.origin) + np.array(self.mre_labels.shape),
        ]

        print(f"MRE image min/maxes: {min_max_img}")

        np.testing.assert_array_less(min_max_img[0], min_max_mesh[0])
        np.testing.assert_array_less(min_max_mesh[1], min_max_img[1])

        # check mesh regions
        shape = (64, 64, 64)
        radii = [5, 10, 15, 20]  # Radii for spheres
        center = np.array(shape) // 2
        origin = (4.0, -1.0, 2.0)

        center_dist_mesh = np.linalg.norm(
            self.mdl.centroid_table - (center + origin), axis=1
        )

        for i in range(1, 5):
            print(i)
            points_in_region = self.mdl.element_table[:, 1] == i

            np.testing.assert_array_less(
                center_dist_mesh[points_in_region], radii[i - 1] + 0.1
            )

    def test_check_mre_inputs(self):
        shape = (64, 64, 64)
        center = np.array(shape) // 2
        radii = [2, 3, 4]
        origin = (4.0, -1.0, 2.0)

        map = spatial_map(self.mre_labels)

        center_dist_mesh = np.linalg.norm(
            map[:, :3] - (center + origin), axis=1
        )

        for i in range(1, 4):
            print(i)
            points_in_region = map[:, 3] == i
            print(center + origin)
            print(map[points_in_region])
            print(map[points_in_region, :3] - (center + origin))
            print(center_dist_mesh[points_in_region])

            np.testing.assert_array_less(
                center_dist_mesh[points_in_region], radii[i - 1] + 0.1
            )

    def test_mapping(self):
        region_props = [[0, 0, 0], [1, 2, 3], [1, 2, 3], [1, 2, 3]]

        test_mdl = map_MRE_to_mesh(
            self.mdl,
            label_img=self.mre_labels,
            target_region_id=1,
            region_properties=region_props,
        )

        # check part ID range
        print(f"unique PIDS: {np.unique(test_mdl.element_table[:, 1])}")
        assert np.max(test_mdl.element_table[:, 1]) == 7

        for i in range(2, 8):
            assert i in test_mdl.element_table[:, 1]

        # check all points are within the regions defined by the MRE labels
        target_region_mask = self.mdl.element_table[:, 1] == 1

        updated_elements = test_mdl.element_table[target_region_mask, :]

        updated_centroids = test_mdl.centroid_table[target_region_mask, :]

        center_coord = self.obj_center + self.obj_origin

        dist = np.linalg.norm(updated_centroids - center_coord, axis=0)

        for idx, r in enumerate(self.mapping_radii):
            submask = updated_elements[:, 1] == idx + 5

            np.testing.assert_array_less(dist[submask], r)

        for point in updated_centroids:
            assert (
                np.linalg.norm(point - self.obj_center + self.obj_origin)
                <= 5.0
            )

        # check material assignments
        assert len(test_mdl.material_info) == 3
