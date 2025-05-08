import pytest

import numpy as np

from MRI2FE import COM_align, point_cloud_spacing


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

        # no data input
        with pytest.raises(ValueError):
            spacing = point_cloud_spacing(dims=dim)

        # wrong dimensions
        with pytest.raises(ValueError):
            spacing = point_cloud_spacing(dims=dim, lims=lim)

        with pytest.raises(ValueError):
            spacing = point_cloud_spacing(dims=dim, points=img)
