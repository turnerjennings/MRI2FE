import pytest

import numpy as np

import os

import meshio

from MRI2FE import mesh_from_nifti

class TestMeshing:

    def verify_radius(self):
        
        root_path = os.getcwd()

        test_file = os.path.join(root_path, 
                                 "test", 
                                 "test_data", 
                                 "test_concentric_spheres.nii")
        
        mesh = mesh_from_nifti(test_file)

        points = mesh.points

        distance = np.linalg.norm(points,axis=1)

        for i in range(distance.shape[0]):
            assert (distance[i] < 20.0)