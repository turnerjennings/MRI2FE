import nibabel as nb
import os

path = "./brainwebmri/"
outpath = "C:/Users/jennings.t/Downloads/brainwebnii/"


for file in os.listdir(path):
    if file.endswith(".mnc"):
        print(f"Processing file {file}")
        bw_mnc = nb.load(os.path.join(path, file))
        bw_model = nb.Nifti1Image(bw_mnc.dataobj, bw_mnc.affine, bw_mnc.header)

        nb.save(bw_model, os.path.join(outpath, file[:-4] + ".nii"))
