import pytest

import numpy as np
import ants

from MRI2FE import (
    COM_align,
    point_cloud_spacing,
    ants_affine,
    spatial_map,
    element_centroids,
)


def generate_test_array():
    arr1 = np.array(  # cube with COM at (0,0,0)
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [
                1.0,
                1.0,
                -1.0,
            ],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [
                1.0,
                1.0,
                1.0,
            ],
        ]
    )

    arr2 = np.array([[2.0, 2.0, 2.0]])

    arr3 = np.array(  # cube with COM at (0,0,0) + random points
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [
                1.0,
                1.0,
                -1.0,
            ],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [
                1.0,
                1.0,
                1.0,
            ],
            [4.0, 5.0, 6.0],
            [2.0, 8.0, 10.0],
        ]
    )

    return arr1, arr2, arr3


class TestCOMAlign:
    def test_com_align(self):
        # create fixed image and moving image

        fi, mi, _ = generate_test_array()

        transformed_arr = COM_align(fixed=fi, moving=mi)

        np.testing.assert_almost_equal(
            transformed_arr, np.array([[0.0, 0.0, 0.0]])
        )

    def test_com_mask(self):
        # create fixed image and moving image

        fi_mask, mi, fi = generate_test_array()

        transformed_arr = COM_align(fixed=fi, moving=mi, fixed_mask=fi_mask)

        np.testing.assert_almost_equal(
            transformed_arr, np.array([[0.0, 0.0, 0.0]])
        )

    def test_raises(self):
        with pytest.raises(ValueError):
            COM_align(fixed=None, moving=np.array([[1, 2, 3]]))
        with pytest.raises(ValueError):
            COM_align(fixed=np.array([[1, 2, 3]]), moving=None)
        with pytest.raises(ValueError):
            COM_align(fixed=np.array([1, 2, 3]), moving=np.array([[1, 2, 3]]))
        with pytest.raises(ValueError):
            COM_align(
                fixed=np.array([[1, 2, 3]]), moving=np.array([[[1, 2, 3]]])
            )
        with pytest.raises(ValueError):
            COM_align(fixed=np.array([[1, 2]]), moving=np.array([[1, 2, 3]]))
        with pytest.raises(ValueError):
            COM_align(
                fixed=np.array([[1, 2, 3]]),
                moving=np.array([[1, 2, 3]]),
                fixed_mask=np.array([1, 2, 3]),
            )
        with pytest.raises(ValueError):
            COM_align(
                fixed=np.array([[1, 2, 3]]),
                moving=np.array([[1, 2, 3]]),
                moving_mask=np.array([[1, 2]]),
            )


class TestPointCloudSpacing:
    def test_spacing_cloud(self):
        # check point cloud spacing called with point coordinates

        img, _, _ = generate_test_array()

        dim = (100, 100, 100)

        spacing = point_cloud_spacing(dims=dim, points=img)

        ref = np.array([0.02, 0.02, 0.02])

        np.testing.assert_almost_equal(spacing, ref)

    def test_spacing_limits(self):
        # check point cloud spacing called with dimension limits

        lims = np.array([[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]])

        dim = (100, 100, 100, 100)

        spacing = point_cloud_spacing(dims=dim, lims=lims)

        ref = np.array([0.02, 0.02, 0.02, 0.02])

        np.testing.assert_almost_equal(spacing, ref)

    def test_raises(self):
        dim = (100, 100, 100, 100)
        lim = np.array([[-1, 1]])
        img, _, _ = generate_test_array()

        with pytest.raises(ValueError):
            point_cloud_spacing(dims=None, points=img)
        with pytest.raises(ValueError):
            point_cloud_spacing(dims=dim)
        with pytest.raises(ValueError):
            point_cloud_spacing(dims=dim, lims=lim)
        with pytest.raises(ValueError):
            point_cloud_spacing(dims=dim, points=img)


class TestAntsAffine:
    def test_raises(self):
        with pytest.raises(ValueError):
            ants_affine(img=None)
        with pytest.raises(TypeError):
            ants_affine(img=np.array([1, 2, 3]))

    def test_ants_affine(self):
        test_data = np.random.rand(10, 10, 10)
        test_img = ants.from_numpy(
            data=test_data,
            origin=(0, 0, 0),
            spacing=(1, 1, 1),
            direction=np.eye(3),
        )

        affine = ants_affine(test_img)
        assert affine.shape == (4, 4)
        np.testing.assert_array_equal(affine[3], [0, 0, 0, 1])


class TestSpatialMap:
    def test_raises(self):
        with pytest.raises(ValueError):
            spatial_map(infile=None)
        with pytest.raises(TypeError):
            spatial_map(infile=np.array([1, 2, 3]))

    def test_spatial_map(self):
        test_data = np.ones((2, 2, 2))
        test_img = ants.from_numpy(
            data=test_data,
            origin=(0, 0, 0),
            spacing=(1, 1, 1),
            direction=np.eye(3),
        )

        result = spatial_map(test_img)
        assert result.shape[1] == 4
        assert len(result) == np.prod(test_data.shape)


class TestElementCentroids:
    def test_centroids(self):
        element = np.array([1, 1, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5])
        coords = np.array(
            [
                [1, 0.0, 0.0, 0.0],
                [2, 1.0, 0.0, 0.0],
                [3, 2.0, 0.0, 0.0],
                [4, 3.0, 0.0, 0.0],
                [5, 4.0, 0.0, 0.0],
            ]
        )

        cx = element_centroids(element, node_coords=coords)

        print(cx)

        np.testing.assert_almost_equal(cx, (2.0, 0.0, 0.0))
