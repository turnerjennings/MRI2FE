import pytest

from ants import from_numpy

import numpy as np

from MRI2FE.Meshing import find_surface

def generate_MRI_sphere():

    array_shape = (101,101,101)

    radius = 20

    spacing = (1.,1.,1.)

    origin = (0.,0.,0.)

    direction = np.eye(3)

    arr = np.zeros(array_shape)


    center = np.array([50, 50, 50])

    # Create coordinate grids
    z, y, x = np.ogrid[:array_shape[0], :array_shape[1], :array_shape[2]]
    
    # Calculate squared distance from center for each point
    # (Using squared distance avoids computing square roots)
    dist_squared = (z - center[0])**2 + (y - center[1])**2 + (x - center[2])**2
    
    # Set points within the radius to 1
    mask = dist_squared <= radius**2
    arr[mask] = 1

    sphere=from_numpy(arr,origin=origin, spacing=spacing,direction=direction)

    return sphere

class TestNiftiToMesh:

    def test_sphere_surface(self):
        image = generate_MRI_sphere()

        verts, faces = find_surface(img=image, low_value=0., high_value=1.)

        print(verts)

        for i in range(verts.shape[0]):
            dist = np.linalg.norm(verts[i,:]-np.array([50.,50.,50.]))
            assert (dist - 20)**2 <= 0.5



