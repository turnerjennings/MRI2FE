from ants.core.ants_image import ANTsImage

from skimage.measure import marching_cubes

def find_surface(img:ANTsImage, low_value:float, high_value:float):

    img_data = img.numpy()

    img_spacing = img.spacing

    if len(img_data.shape) > 3:
        raise ValueError("Image must be 3D")
    
    threshold = 0.5*(high_value - low_value) + low_value
    
    verts, faces, _, _ = marching_cubes(volume=img_data,
                                        level=threshold,
                                        spacing=img_spacing)

    return verts, faces
    
