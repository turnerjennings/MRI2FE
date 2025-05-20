#define CATCH_CONFIG_MAIN
#include<catch2/catch_test_macros.hpp>
#include<catch2/catch_approx.hpp>
#include "meshing.hpp"

/**
 * Test basic add function to confirm CMAKE format and compilation
 */
TEST_CASE("test_debug_add") {
    REQUIRE(_debug_add(1,2) == 3);
    REQUIRE(_debug_add(0,0) == 0);
    REQUIRE(_debug_add(-1,2) == 1);
}