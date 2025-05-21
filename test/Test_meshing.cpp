#define CATCH_CONFIG_MAIN
#include<catch2/catch_test_macros.hpp>
#include<catch2/catch_approx.hpp>
#include "meshing.hpp"
#include "io.hpp"
#include<vtkImageData.h>

#include <string>
#include <filesystem>

namespace fs=std::filesystem;

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
    std::string fname = "zstat1.nii";

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
    REQUIRE(dimensions[2] == 21);

    //check spacing
    double spacing[3];
    test_data->GetSpacing(spacing);

    REQUIRE(spacing[0] == Catch::Approx(4.0));
    REQUIRE(spacing[1] == Catch::Approx(4.0));
    REQUIRE(spacing[2] == Catch::Approx(6.0));


}

TEST_CASE("load nifti throws"){
    std::string resource_dir = TEST_RESOURCE_DIR;
    std::string fname = "zstat1.gz";

    fs::path full_path = fs::path(resource_dir) / fname;
    std::string fpath = full_path.string();

    CAPTURE(fpath);

    REQUIRE_THROWS(load_nifti(fpath));


    fname = "zstat2.nii";

    full_path = fs::path(resource_dir) / fname;
    fpath = full_path.string();

    REQUIRE_THROWS(load_nifti(fpath));
}