#define CATCH_CONFIG_MAIN
#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "meshing_core.hpp"

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include<CGAL/Mesh_triangulation_3.h>
#include<CGAL/Mesh_complex_3_in_triangulation_3.h>
#include<CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/IO/File_binary_mesh_3.h>

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


TEST_CASE("Meshing no lloyd"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.inr";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    Mesh_criteria criteria(params::facet_angle(30).facet_size(1).facet_distance(4).
    cell_radius_edge_ratio(3).cell_size(1));

    std::string outpath = mesh_model(fpath, criteria, false);

    INFO(outpath);

    fs::path check_path = outpath;

    REQUIRE(fs::exists(check_path));
}

TEST_CASE("Meshing with lloyd"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.inr";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    Mesh_criteria criteria(params::facet_angle(30).facet_size(1).facet_distance(4).
    cell_radius_edge_ratio(3).cell_size(1));

    std::string outpath = mesh_model(fpath, criteria, true);

    INFO(outpath);

    fs::path check_path = outpath;

    REQUIRE(fs::exists(check_path));
}

TEST_CASE("Mesh_wrapper"){

    //load nifti
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "test_concentric_spheres.inr";

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