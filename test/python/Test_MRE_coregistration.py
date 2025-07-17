import pytest
import numpy as np

from ants import get_ants_data, image_read, registration

from MRI2FE.MRE import coregister_MRE_images


class TestMRECoreg:
    def test_gstar_from_buffer(self):
        geom = image_read(get_ants_data("r16"))

        gp = image_read(get_ants_data("r27"))
        gpp = image_read(get_ants_data("r27"))

        test_dict = coregister_MRE_images(segmented_geom=geom, gp_list=gp, gpp_list=gpp)

        assert "gp" in test_dict
        assert "gpp" in test_dict
        assert "transform" in test_dict

        print(test_dict["gp"])
        print(test_dict["gpp"])

        print(test_dict["gp"].min(), test_dict["gp"].max())

        np.testing.assert_almost_equal(
            test_dict["gp"].numpy(), test_dict["gpp"].numpy(), decimal=2
        )

    def test_gstar_from_buffer_list(self):
        geom = image_read(get_ants_data("r16"))

        gp = image_read(get_ants_data("r27"))
        gpp = image_read(get_ants_data("r27"))

        test_dict = coregister_MRE_images(
            segmented_geom=geom, gp_list=[gp, gp, gp], gpp_list=[gpp, gpp, gpp]
        )

        for out in test_dict:
            assert "gp" in out
            assert "gpp" in out
            assert "transform" in out

            print(out["gp"])
            print(out["gpp"])

            print(out["gp"].min(), out["gp"].max())

            np.testing.assert_almost_equal(
                out["gp"].numpy(), out["gpp"].numpy(), decimal=2
            )

    def test_ddsr_from_fpath(self):
        geom = get_ants_data("r16")

        mu = get_ants_data("r27")
        xi = get_ants_data("r27")

        test_dict = coregister_MRE_images(
            segmented_geom=geom, mu_list=[mu], xi_list=[xi]
        )

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
            test_dict = coregister_MRE_images(
                segmented_geom=geom, gp_list=mu, xi_list=xi
            )

        with pytest.raises(ValueError):
            test_dict = coregister_MRE_images(
                segmented_geom=geom, mu_list=mu, gpp_list=xi
            )

        with pytest.raises(FileNotFoundError):
            test_dict = coregister_MRE_images(segmented_geom="nonexistent_file.nii.gz")

        with pytest.raises(FileNotFoundError):
            test_dict = coregister_MRE_images(
                segmented_geom=geom,
                geom_mask="nonexistent_mask.nii.gz",
                gp_list=mu,
                gpp_list=xi,
            )

        with pytest.raises(TypeError):
            test_dict = coregister_MRE_images(segmented_geom=123)

        with pytest.raises(TypeError):
            test_dict = coregister_MRE_images(segmented_geom=geom, geom_mask=[1, 2, 3])

        with pytest.raises(ValueError):
            test_dict = coregister_MRE_images(
                segmented_geom=geom,
                gp_list=mu,
                gpp_list=xi,
                imgout="/nonexistent_directory/output",
            )
