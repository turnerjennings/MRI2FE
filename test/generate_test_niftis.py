import numpy as np
from ants import from_numpy, image_write
import os

# Parameters
shape = (64, 64, 64)
center = np.array(shape) // 2
radii = [5, 10, 15, 20]  # Radii for spheres
labels = list(range(1, len(radii)+1))  # Labels: 1, 2, 3, ...

spacing=(1.,1.,1.)
origin=(4.,-1.,2.)

# Create empty volume
volume = np.zeros(shape, dtype=np.uint8)

# Grid of coordinates
zz, yy, xx = np.meshgrid(
    np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing='ij'
)
dist = np.sqrt((xx - center[2])**2 + (yy - center[1])**2 + (zz - center[0])**2)

# Assign labels to voxels inside each sphere
for i, r in enumerate(radii):
    if i == 0:
        mask = dist <= r
    else:
        mask = (dist <= r) & (dist > radii[i-1])
    volume[mask] = labels[i]

# Save as NIfTI
nifti_img = from_numpy(data=volume, 
                       origin=origin,
                       spacing=spacing)

cwd = os.getcwd()

output_dir = os.path.join(cwd,"test","test_data")
os.makedirs(output_dir,exist_ok=True)

image_write(nifti_img, os.path.join(output_dir,"test_concentric_spheres.nii"))