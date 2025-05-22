#define CATCH_CONFIG_MAIN
#include<catch2/catch_test_macros.hpp>
#include<catch2/catch_approx.hpp>
#include "meshing.hpp"
#include "io.hpp"
#include<vtkImageData.h>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include<CGAL/Mesh_triangulation_3.h>
#include<CGAL/Mesh_complex_3_in_triangulation_3.h>
#include<CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/IO/File_binary_mesh_3.h>


#include<CGAL/IO/read_vtk_image_data.h>

#include <string>
#include <filesystem>

namespace fs=std::filesystem;

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

/**
 * Test basic add function to confirm CMAKE format and compilation
 */
TEST_CASE("test_debug_add") {
    REQUIRE(_debug_add(1,2) == 3);
    REQUIRE(_debug_add(0,0) == 0);
    REQUIRE(_debug_add(-1,2) == 1);
}

TEST_CASE("load nifti") {
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.nii";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    vtkImageData* test_data = load_nifti(fpath);

    //check dimensions
    int dimensions[3];
    test_data->GetDimensions(dimensions);

    REQUIRE(test_data != nullptr);

    REQUIRE(dimensions[0] == 64);
    REQUIRE(dimensions[1] == 64);
    REQUIRE(dimensions[2] == 64);

    //check spacing
    double spacing[3];
    test_data->GetSpacing(spacing);

    REQUIRE(spacing[0] == Catch::Approx(1.0));
    REQUIRE(spacing[1] == Catch::Approx(1.0));
    REQUIRE(spacing[2] == Catch::Approx(1.0));

    //check origin
    double origin[3];
    test_data->GetOrigin(origin);

    REQUIRE(origin[0] == Catch::Approx(-4.0));
    REQUIRE(origin[1] == Catch::Approx(1.0));
    REQUIRE(origin[2] == Catch::Approx(2.0));   


}

TEST_CASE("load nifti throws"){
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "not_real.txt";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    REQUIRE_THROWS(load_nifti(fpath));


    fname = "not_real.nii";

    full_path = fs::path(resource_dir) / fname;
    fpath = full_path.string();

    REQUIRE_THROWS(load_nifti(fpath));
}

TEST_CASE("Meshing no lloyd"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.nii";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    vtkImageData* test_data = load_nifti(fpath);

    Mesh_criteria criteria(params::facet_angle(30).facet_size(1).facet_distance(4).
    cell_radius_edge_ratio(3).cell_size(1));

    std::string outpath = mesh_model(test_data, criteria, false);

    INFO(outpath);

    fs::path check_path = outpath;

    REQUIRE(fs::exists(check_path));
}

TEST_CASE("Meshing with lloyd"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.nii";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    vtkImageData* test_data = load_nifti(fpath);

    Mesh_criteria criteria(params::facet_angle(30).facet_size(1).facet_distance(4).
    cell_radius_edge_ratio(3).cell_size(1));

    std::string outpath = mesh_model(test_data, criteria, true);

    INFO(outpath);

    fs::path check_path = outpath;

    REQUIRE(fs::exists(check_path));
}

TEST_CASE("Mesh_wrapper"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.nii";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);


    std::string outpath = mesh_wrapper(fpath,
    false,
    30.0,
    1.0,
    4.0,
    3.0,
    1.0);

    INFO(outpath);

    fs::path check_path = outpath;

    REQUIRE(fs::exists(check_path));
}