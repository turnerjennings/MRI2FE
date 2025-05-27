import pytest
import numpy as np

from ants import get_ants_data, image_read, registration

from MRI2FE.MRE import coregister_MRE_images


class TestMRECoreg:
    def test_gstar_from_buffer(self):
        geom = image_read(get_ants_data("r16"))

        gp = image_read(get_ants_data("r27"))
        gpp = image_read(get_ants_data("r27"))

        test_dict = coregister_MRE_images(geom=geom, gp=gp, gpp=gpp)

        assert "gp" in test_dict
        assert "gpp" in test_dict
        assert "transform" in test_dict

        print(test_dict["gp"])
        print(test_dict["gpp"])

        print(test_dict["gp"].min(), test_dict["gp"].max())

        np.testing.assert_almost_equal(
            test_dict["gp"].numpy(), test_dict["gpp"].numpy(), decimal=2
        )

    def test_ddsr_from_fpath(self):
        geom = get_ants_data("r16")

        mu = get_ants_data("r27")
        xi = get_ants_data("r27")

        test_dict = coregister_MRE_images(geom=geom, mu=mu, xi=xi)

        assert "mu" in test_dict
        assert "xi" in test_dict
        assert "transform" in test_dict

        print(test_dict["mu"])
        print(test_dict["xi"])

        print(test_dict["mu"].min(), test_dict["mu"].max())

        np.testing.assert_almost_equal(
            test_dict["mu"].numpy(), test_dict["xi"].numpy(), decimal=2
        )

    def test_raises(self):
        geom = get_ants_data("r16")

        mu = get_ants_data("r27")
        xi = get_ants_data("r27")

        with pytest.raises(ValueError):
            test_dict = coregister_MRE_images(geom=geom, gp=mu, xi=xi)

            test_dict = coregister_MRE_images(geom=geom, mu=mu, gpp=xi)

        with pytest.raises(TypeError):
            test_dict = coregister_MRE_images(geom_mask=geom)
