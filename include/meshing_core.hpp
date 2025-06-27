#pragma once

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include<CGAL/Mesh_triangulation_3.h>
#include<CGAL/Mesh_complex_3_in_triangulation_3.h>
#include<CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include<filesystem>
#include<string>

namespace fs=std::filesystem;
namespace params = CGAL::parameters;

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh_domain = CGAL::Labeled_mesh_domain_3<K>;
 
// Triangulation
using Tr   = CGAL::Mesh_triangulation_3<Mesh_domain,
                                        CGAL::Default,
                                        CGAL::Parallel_if_available_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
 
// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

int _debug_add(int i, int j);

std::string mesh_model(std::string fpath, 
    const Mesh_criteria criteria, 
    const bool lloyd
    );

std::string mesh_wrapper(std::string filePath,
                        const bool optimize,
                        const double facetAngle,
                        const double facetSize,
                        const double facetDistance,
                        const double cellRadiusEdgeRatio,
                        const double cellSize); 