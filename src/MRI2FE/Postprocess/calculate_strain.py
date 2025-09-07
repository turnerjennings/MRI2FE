import numpy as np
from ants import from_numpy
from ants.core.ants_image import ANTsImage

from ..utilities import check_xyz


def MRI_strain(
    img_1: ANTsImage,
    img_2: ANTsImage,
    img_3: ANTsImage,
) -> ANTsImage:
    """Calculate the strain field for an MRI object given displacement in three principal directions

    Args:
        img_1 (ANTsImage): Displacement field in first principal direction, shape (l,m,n)
        img_2 (ANTsImage): Displacement field in second principal direction, shape (l,m,n)
        img_3 (ANTsImage): Displacement field in third principal direction, shape (l,m,n)

    Raises:
        TypeError: If inputs are not ANTsImage objects
        ValueError: If inputs are None

    Returns:
        AntsImage: Strain field image, shape (l,m,n,6) with strain ordered e_11, e_22, e_33, e_12, e_13, e_23
    """
    if img_1 is None:
        raise ValueError("img_1 cannot be None")
    if img_2 is None:
        raise ValueError("img_2 cannot be None")
    if img_3 is None:
        raise ValueError("img_3 cannot be None")

    if not isinstance(img_1, ANTsImage):
        raise TypeError("img_1 must be an ANTsImage object")
    if not isinstance(img_2, ANTsImage):
        raise TypeError("img_2 must be an ANTsImage object")
    if not isinstance(img_3, ANTsImage):
        raise TypeError("img_3 must be an ANTsImage object")

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

    strain_img = from_numpy(
        data=e_out,
        spacing=img_1.spacing,
        direction=img_1.direction,
        origin=img_1.origin,
        has_components=True,
    )

    return strain_img
