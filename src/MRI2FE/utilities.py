import numpy as np
import ants
from typing import Union


def COM_align(
    fixed: np.ndarray,
    moving: np.ndarray,
    fixed_mask: np.ndarray = None,
    moving_mask: np.ndarray = None,
) -> np.ndarray:
    """Align the centers of mass of two point clouds, operating either on the entire point cloud or on a masked region.

    Args:
        fixed (np.ndarray): Reference point cloud to be aligned to.
        moving (np.ndarray): Point cloud to align with fixed cloud.
        fixed_mask (np.ndarray, optional): Fixed point cloud subset to define COM by. Defaults to None.
        moving_mask (np.ndarray, optional): Moving point cloud subset to define COM by. Defaults to None.

    Raises:
        TypeError: If inputs are not numpy arrays when provided
        ValueError: If required inputs are None or have invalid dimensions

    Returns:
        np.ndarray: moving point cloud rigidly aligned to the fixed point cloud COM
    """
    # Validate required inputs are not None
    if fixed is None:
        raise ValueError("fixed cannot be None")
    if moving is None:
        raise ValueError("moving cannot be None")

    # confirm numpy array input and validate dimensions
    fixed = np.asarray(fixed)
    if len(fixed.shape) != 2:
        raise ValueError("fixed must be a 2D array")

    moving = np.asarray(moving)
    if len(moving.shape) != 2:
        raise ValueError("moving must be a 2D array")

    # Validate matching dimensions
    if fixed.shape[1] != moving.shape[1]:
        raise ValueError(
            "fixed and moving must have the same number of dimensions"
        )

    # Validate masks if provided
    if fixed_mask is not None:
        fixed_mask = np.asarray(fixed_mask)
        if len(fixed_mask.shape) != 2:
            raise ValueError("fixed_mask must be a 2D array")
        if fixed_mask.shape[1] != fixed.shape[1]:
            raise ValueError(
                "fixed_mask must have same number of dimensions as fixed"
            )

    if moving_mask is not None:
        moving_mask = np.asarray(moving_mask)
        if len(moving_mask.shape) != 2:
            raise ValueError("moving_mask must be a 2D array")
        if moving_mask.shape[1] != moving.shape[1]:
            raise ValueError(
                "moving_mask must have same number of dimensions as moving"
            )

    # Calculate Centers of Mass
    if fixed_mask is not None:
        COM_fixed = np.mean(fixed_mask, axis=0)
    else:
        COM_fixed = np.mean(fixed, axis=0)

    if moving_mask is not None:
        COM_moving = np.mean(moving_mask, axis=0)
    else:
        COM_moving = np.mean(moving, axis=0)

    # calculate translation and create transformation matrix
    offset = COM_fixed - COM_moving

    transform = np.array(
        [
            [1, 0, 0, offset[0]],
            [0, 1, 0, offset[1]],
            [0, 0, 1, offset[2]],
            [0, 0, 0, 1],
        ]
    )

    # apply transform
    moving_augment = np.hstack((moving, np.ones((moving.shape[0], 1))))

    moving_transformed = moving_augment @ transform.T

    return moving_transformed[:, 0:3]


def point_cloud_spacing(
    dims: Union[list, tuple, np.ndarray],
    points: Union[list, tuple, np.ndarray] = None,
    lims: Union[list, tuple, np.ndarray] = None,
) -> tuple:
    """Return the voxel spacing necessary to cover a point cloud with given voxel dimensions
    Can be called by either providing a point cloud or directly providing field limits.

    Args:
        dims (tuple): tuple of voxel count along each of m dimensions
        points (Union[list, tuple, np.ndarray], optional): locations of n points in m dimensions with shape (n,m)
        lims (np.ndarray, optional): Min/max of point cloud space in m dimensions with shape (m,2)

    Raises:
        TypeError: If inputs are not of correct type
        ValueError: If dims is None or if neither points nor lims is provided

    Returns:
        tuple: spacing in each dimension
    """
    # Validate dims is not None
    if dims is None:
        raise ValueError("dims cannot be None")

    # check that either points or lims is provided
    if points is None and lims is None:
        raise ValueError("Must provide either points or lims for calculation")

    # convert inputs to np.ndarray
    dims = np.asarray(dims)

    if points is not None:
        points = np.asarray(points)

    if lims is not None:
        lims = np.asarray(lims)

    # check that dimensions match
    if points is not None and not (dims.shape[0] == points.shape[-1]):
        raise ValueError(
            "dims must have the same dimension as the point cloud"
        )

    if lims is not None and not (dims.shape[0] == lims.shape[0]):
        raise ValueError("dims must have the same dimension as lims")

    # find limits of point cloud if not provided
    if lims is None:
        mins = np.min(points, axis=0)
        maxes = np.max(points, axis=0)
    else:
        mins = lims[:, 0]
        maxes = lims[:, 1]

    delta = maxes - mins

    spacing = delta / dims

    return spacing


def ants_affine(img: ants.core.ants_image.ANTsImage) -> np.ndarray:
    """Extract affine transformation matrix from ANTsImage.

    Args:
        img (ants.core.ants_image.ANTsImage): Input ANTs image

    Raises:
        TypeError: If img is not an ANTsImage
        ValueError: If img is None

    Returns:
        np.ndarray: 4x4 affine transformation matrix
    """
    # Validate input is not None
    if img is None:
        raise ValueError("img cannot be None")

    # Validate input type
    if not isinstance(img, ants.core.ants_image.ANTsImage):
        raise TypeError("img must be an ANTsImage object")

    # Extract image properties
    spacing = np.array(img.spacing)  # (sx, sy, sz)
    direction = np.array(img.direction).reshape((3, 3))  # 3x3 rotation matrix
    origin = np.array(img.origin)  # (ox, oy, oz)

    # Compute affine matrix
    affine_matrix = np.eye(4)  # Initialize as identity
    affine_matrix[:3, :3] = direction * spacing  # Scale direction by spacing
    affine_matrix[:3, 3] = origin  # Set translation

    return affine_matrix


def check_xyz(
    a: ants.core.ants_image.ANTsImage,
    b: ants.core.ants_image.ANTsImage,
    c: ants.core.ants_image.ANTsImage,
) -> bool:
    # Check if any of the shapes differ
    if not (a.shape == b.shape and b.shape == c.shape):
        return True

    # Check if any directions differ
    # Use np.array_equal for comparing arrays
    if not (
        np.array_equal(a.direction, b.direction)
        and np.array_equal(b.direction, c.direction)
    ):
        return True

    # Check if any origins differ
    # Use np.array_equal for comparing arrays
    if not (
        np.array_equal(a.origin, b.origin)
        and np.array_equal(b.origin, c.origin)
    ):
        return True

    # Check if any spacings differ
    # Convert to tuples if they're numpy arrays
    a_spacing = (
        tuple(a.spacing) if isinstance(a.spacing, np.ndarray) else a.spacing
    )
    b_spacing = (
        tuple(b.spacing) if isinstance(b.spacing, np.ndarray) else b.spacing
    )
    c_spacing = (
        tuple(c.spacing) if isinstance(c.spacing, np.ndarray) else c.spacing
    )

    if not (a_spacing == b_spacing and b_spacing == c_spacing):
        return True

    # If we get here, none of the conditions for being different were met
    return False


def spatial_map(infile: ants.core.ants_image.ANTsImage) -> np.ndarray:
    """Convert data from NIFTI file to (4,n) array of x,y,z,voxel value in physical space

    Args:
        infile (nb.nifti1): NIFTI file to convert

    Raises:
        TypeError: If infile is not an ANTsImage
        ValueError: If infile is None

    Returns:
        coordinates_and_values (np.ndarray): array of x,y,z coordinates in physical space with the corresponding voxel value
    """
    # Validate input is not None
    if infile is None:
        raise ValueError("infile cannot be None")

    # Validate input type
    if not isinstance(infile, ants.core.ants_image.ANTsImage):
        raise TypeError("infile must be an ANTsImage object")

    # extract data from NIFTI
    image_data = infile.numpy()
    affine = ants_affine(infile)
    dimensions = image_data.shape

    # create coordinates in voxel space
    x = np.arange(dimensions[0])
    y = np.arange(dimensions[1])
    z = np.arange(dimensions[2])
    xv, yv, zv = np.meshgrid(x, y, z, indexing="ij")

    voxel_coords = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
    voxel_coords_homogeneous = np.hstack(
        [voxel_coords, np.ones((voxel_coords.shape[0], 1))]
    )

    # map voxel coordinates to physical space
    physical_coords = voxel_coords_homogeneous @ affine.T
    voxel_values = np.round(image_data.ravel())
    coordinates_and_values = np.hstack(
        [physical_coords[:, :3], voxel_values[:, np.newaxis]]
    )

    return coordinates_and_values


def element_centroids(elnodes:Union[np.ndarray,tuple,list], 
                      node_coords:np.ndarray) -> np.ndarray:
    """Calculate the centroid of an element

    Args:
        elnodes (np.array): 1D array containing EID, PID, and connected nodes for arbitrary 10-node element
        node_coords (np.array): 2D array containing NID, xyz coordinates

    Raises:
        TypeError: If inputs are not numpy arrays
        ValueError: If array dimensions or contents are invalid

    Returns:
        centroid (np.array): (3,) array containing average coordinate of the element
    """
    # Validate input types
    if not isinstance(elnodes, (np.ndarray,tuple,list)):
        raise TypeError("elnodes must be a numpy array, tuple, or list")
    if not isinstance(node_coords, np.ndarray):
        raise TypeError("node_coords must be a numpy array")

    # Validate array dimensions
    if len(elnodes.shape) != 1:
        raise ValueError("elnodes must be a 1D array")
    if len(node_coords.shape) != 2:
        raise ValueError("node_coords must be a 2D array")

    # Validate array sizes
    if elnodes.shape[0] < 12:  # Must have at least EID, PID, and 10 nodes
        raise ValueError(
            "elnodes must contain at least 12 elements (EID, PID, and 10 nodes)"
        )
    if node_coords.shape[1] != 4:  # Must have NID and xyz coordinates
        raise ValueError("node_coords must have 4 columns (NID, x, y, z)")

    nodes = np.unique(elnodes[2:])

    coords=[]
    for node in nodes:
        if node not in node_coords[:,0]:
            raise ValueError(f"Point {node} not found in node array during centroid calculation")
        coords.append(node_coords[node_coords[:,0] == node,1:])
    coords = np.array(coords)

    n = len(coords)
    if n < 2:
        raise ValueError(f"only {n} points found during centroid calculation")
    
    cx = np.mean(coords,axis=0).flatten()

    return cx
        
