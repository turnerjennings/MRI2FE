import numpy as np
import meshio
import sys
import os

# Add the src directory to the path to use the local version
sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "src")
)

from MRI2FE.models.femodel import FEModel


# Create a mock element_centroids function to avoid import issues
def mock_element_centroids(*args, **kwargs):
    return np.array([0, 0, 0])


# Create a mock utilities module
class MockUtilities:
    element_centroids = mock_element_centroids


# Add the mock to sys.modules
sys.modules["MRI2FE.utilities"] = MockUtilities()
sys.modules["MRI2FE.utilities.element_centroids"] = MockUtilities()


def test_from_meshio_appends_data():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    cells = [("tetra", np.array([[0, 1, 2, 3]]))]
    cell_data = {"medit:ref": [np.array([1])]}
    mesh = meshio.Mesh(points, cells, cell_data=cell_data)

    femodel = FEModel()
    femodel.from_meshio(mesh, element_type="tetra", region_names=["part1"])
    femodel.from_meshio(mesh, element_type="tetra", region_names=["part2"])

    assert femodel.node_table.shape[0] == 8
    assert femodel.element_table.shape[0] == 2
    assert femodel.metadata["num_nodes"] == 8
    assert femodel.metadata["num_elements"] == 2
    assert "part1" in [p["name"] for p in femodel.part_info.values()]
    assert "part2" in [p["name"] for p in femodel.part_info.values()]
