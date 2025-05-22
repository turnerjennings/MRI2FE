#include "io.hpp"

namespace fs = std::filesystem;

vtkImageData* load_nifti(std::string filepath){

    vtkImageData* image = nullptr;

    //create path object from string and check file type
    fs::path img(filepath);

    if (!fs::exists(img)) {
        throw std::invalid_argument("Input filepath does not exist");
    }

    if (img.extension() != ".nii"){
        throw std::invalid_argument("Input file is not NIFTI");
    }

    //initialize reader
    vtkNIFTIImageReader* reader = vtkNIFTIImageReader::New();
    
    std::string path_string = img.string();
    reader->SetFileName(path_string.c_str());
    reader->Update();

    image = reader->GetOutput();

    // set matrix
    vtkMatrix4x4* mat = reader->GetQFormMatrix();
    if (!mat) {
        mat = reader->GetSFormMatrix();
    }

    // set origin
    if (mat) {
        double origin[3] = {
            mat->GetElement(0, 3),
            mat->GetElement(1, 3),
            mat->GetElement(2, 3)
        };
        image->SetOrigin(origin);
}

    return image;
}