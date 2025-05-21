#pragma once

#include<CGAL/Mesh_triangulation_3.h>
#include<CGAL/Mesh_criteria_3.h>
#include<CGAL/Image_3.h>
#include<CGAL/IO/read_vtk_image_data.h>


#include<pybind11/pybind11.h>

#include<filesystem>

#include<string>

namespace fs=std::filesystem;

int _debug_add(int i, int j);

std::string mesh_model(const vtkImageData* image);