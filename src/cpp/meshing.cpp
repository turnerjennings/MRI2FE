#include "meshing.hpp"

#include<pybind11/pybind11.h>

namespace py=pybind11;

int _debug_add(int i, int j){
    return i + j;
}

std::string mesh_model(std::string fpath, 
    const Mesh_criteria criteria, 
    const bool lloyd
    ){
    
    CGAL::Image_3 image;

    if(!image.read(fpath)){
        throw std::exception("Cannot read file");
    }

    //define domain
    Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image);
    //mesh

    C3t3 c3t3;

    if(lloyd == false){
        c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    } else {
        c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            params::lloyd(params::time_limit(30)).
            no_perturb().
            exude(params::time_limit(10).sliver_bound(10)));
    }

    //write to temp dir
    fs::path tempDir = fs::temp_directory_path();

    fs::path tempFile = tempDir / "tempmesh.mesh";

    std::ofstream outFile(tempFile);
    if (!outFile) {
        throw std::runtime_error("Failed to create temporary file");
    }

    CGAL::IO::write_MEDIT(outFile, c3t3);
    outFile.close();

    return tempFile.string();
}

std::string mesh_wrapper(std::string filePath,
const bool optimize,
const double facetAngle,
const double facetSize,
const double facetDistance,
const double cellRadiusEdgeRatio,
const double cellSize){

    Mesh_criteria criteria(params::facet_angle(facetAngle).
    facet_size(facetSize).
    facet_distance(facetDistance).
    cell_radius_edge_ratio(cellRadiusEdgeRatio).
    cell_size(cellSize));

    std::string outpath = mesh_model(filePath, criteria, optimize);

    if (!fs::exists(fs::path(outpath))) {
        throw std::invalid_argument("Output file not found");
    }

    return outpath;
}

PYBIND11_MODULE(_MESHUTILS, m) {
    m.doc() = "meshing binding for tetrahedral meshing of segmented MRI";

    m.def("mesh_wrapper", 
        &mesh_wrapper,
        py::arg("filePath"),
        py::arg("optimize"),
        py::arg("facetAngle"),
        py::arg("facetSize"),
        py::arg("facetDistance"),
        py::arg("cellRadiusEdgeRatio"),
        py::arg("cellSize"),
    "create a mesh from a segmented NIFTI image.");
}