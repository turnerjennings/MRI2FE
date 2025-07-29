import pytest
from MRI2FE.Postprocess import MRI_strain

import ants
import numpy as np


def image_setup():
    # generate test data
    # assuming u1=xy, u2=yz, u3=xz

    # define shape and spacing
    shp = (100, 100, 100)
    spac = (1.0, 1.0, 1.0)
    dr = np.eye(3)
    orgn = (0.0, 0.0, 0.0)

    i, j, k = np.indices(shp)

    u1 = i * j
    u2 = j * k
    u3 = i * k
    print(f"data shape: {u1.shape}")

    # define output comparison
    e11 = j
    e22 = k
    e33 = i
    e12 = 0.5 * i
    e13 = 0.5 * k
    e23 = 0.5 * j

    e_compare = np.stack([e11, e22, e33, e12, e13, e23], axis=-1)

    # create displacement images
    img_u1 = ants.from_numpy(data=u1, spacing=spac, direction=dr, origin=orgn)

    img_u2 = ants.from_numpy(data=u2, spacing=spac, direction=dr, origin=orgn)

    img_u3 = ants.from_numpy(data=u3, spacing=spac, direction=dr, origin=orgn)

    return img_u1, img_u2, img_u3, e_compare


class TestMRIStrain:
    def test_MRI_strain(self):
        img_u1, img_u2, img_u3, e_compare = image_setup()

        e_out = MRI_strain(img_u1, img_u2, img_u3)

        np.testing.assert_almost_equal(e_out.numpy(), e_compare)

    def test_raises(self):
        # check inconsistent spacing
        img_u1, img_u2, img_u3, _ = image_setup()

        img_u1.set_spacing((2.0, 1.0, 1.0))

        with pytest.raises(ValueError):
            _ = MRI_strain(img_u1, img_u2, img_u3)

        # check inconsistent origin
        img_u1, img_u2, img_u3, _ = image_setup()

        img_u1.set_origin((-1.0, -1.0, -1.0))

        with pytest.raises(ValueError):
            _ = MRI_strain(img_u1, img_u2, img_u3)

        # check inconsistent direction
        img_u1, img_u2, img_u3, _ = image_setup()

        img_u1.set_direction(np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]))

        with pytest.raises(ValueError):
            _ = MRI_strain(img_u1, img_u2, img_u3)

        with pytest.raises(ValueError):
            _ = MRI_strain(None, img_u2, img_u3)
        with pytest.raises(ValueError):
            _ = MRI_strain(img_u1, None, img_u3)
        with pytest.raises(ValueError):
            _ = MRI_strain(img_u1, img_u2, None)

        invalid_array = np.random.rand(100, 100, 100)
        with pytest.raises(TypeError):
            _ = MRI_strain(invalid_array, img_u2, img_u3)
        with pytest.raises(TypeError):
            _ = MRI_strain(img_u1, invalid_array, img_u3)
        with pytest.raises(TypeError):
            _ = MRI_strain(img_u1, img_u2, invalid_array)
