#pragma once

#include <pybind11/pybind11.h>

#include <vtkNew.h>
#include <vtkNIFTIImageReader.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <string>
#include <filesystem>

vtkImageData* load_nifti(std::string filepath);

