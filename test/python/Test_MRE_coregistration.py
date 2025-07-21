import pytest
import numpy as np

from ants import get_ants_data, image_read, threshold_image

from MRI2FE.MRE import coregister_MRE_images, segment_MRE_regions


class TestMRECoreg:
    @pytest.fixture(autouse=True)
    def problem_setup(self):
        self.geom = image_read(get_ants_data("r16"))
        self.geom_mask = threshold_image(
            self.geom,
            low_thresh=np.min(np.nonzero(self.geom.numpy())),
            high_thresh=np.max(self.geom.numpy()),
        )
        print(self.geom)
        self.MRE_geom = image_read(get_ants_data("r27"))
        self.MRE_mask = threshold_image(
            self.MRE_geom,
            low_thresh=np.min(np.nonzero(self.MRE_geom.numpy())),
            high_thresh=np.max(self.geom.numpy()),
        )
        print(self.MRE_geom)
        self.gp = image_read(get_ants_data("r27"))
        self.gpp = image_read(get_ants_data("r27"))

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

        assert len(test_dict["1"]) == 5
        assert len(test_dict["2"]) == 5

        for list in test_dict["1"]:
            assert len(list) == 1

        for list in test_dict["2"]:
            assert len(list) == 1

        assert test_img.max() == 4

    def test_seg_multiple(self):
        ss = image_read("test/test_data/test_stiffness.nii")

        dr = image_read("test/test_data/test_damping_ratio.nii")

        img_list = [(ss, dr), (ss, dr), (ss, dr)]

        test_img, test_dict = segment_MRE_regions(img_list=img_list, n_segs=5)

        assert len(test_dict["1"]) == 5
        assert len(test_dict["2"]) == 5

        for list in test_dict["1"]:
            assert len(list) == 3

        for list in test_dict["2"]:
            assert len(list) == 3

        assert test_img.max() == 4
