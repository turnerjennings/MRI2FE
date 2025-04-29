import ants
import numpy as np
from ..utilities import check_xyz


def MRI_strain(
    img_1: ants.core.ants_image.ANTsImage,
    img_2: ants.core.ants_image.ANTsImage,
    img_3: ants.core.ants_image.ANTsImage,
) -> ants.core.ants_image.ANTsImage:
    # check shape and spacing
    if check_xyz(img_1, img_2, img_3):
        raise ValueError(
            "Images do not share the same spacing, shape, or direction..."
        )

    # extract numpy data
    u1 = img_1.numpy()
    u2 = img_2.numpy()
    u3 = img_3.numpy()

    # calculate gradients
    u11 = np.gradient(u1, img_1.spacing[0], axis=0)
    u12 = np.gradient(u1, img_1.spacing[1], axis=1)
    u13 = np.gradient(u1, img_1.spacing[2], axis=2)

    u21 = np.gradient(u2, img_2.spacing[0], axis=0)
    u22 = np.gradient(u2, img_2.spacing[1], axis=1)
    u23 = np.gradient(u2, img_2.spacing[2], axis=2)

    u31 = np.gradient(u3, img_3.spacing[0], axis=0)
    u32 = np.gradient(u3, img_3.spacing[1], axis=1)
    u33 = np.gradient(u3, img_3.spacing[2], axis=2)

    # calculate strain
    e11 = u11
    e22 = u22
    e33 = u33
    e12 = 0.5 * (u12 + u21)
    e13 = 0.5 * (u13 + u31)
    e23 = 0.5 * (u23 + u32)

    e_out = np.stack([e11, e22, e33, e12, e13, e23], axis=-1)

    strain_img = ants.from_numpy(
        data=e_out,
        spacing=img_1.spacing,
        direction=img_1.direction,
        origin=img_1.origin,
        has_components=True,
    )

    return strain_img
