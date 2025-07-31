import pytest
from MRI2FE import FEModel
import numpy as np
import os


def cube_nodes() -> np.ndarray:
    return np.array(
        [
            [1, 0.0, 0.0, 0.0],
            [2, 1.0, 0.0, 0.0],
            [3, 1.0, 1.0, 0.0],
            [4, 0.0, 1.0, 0.0],
            [5, 0.0, 0.0, 1.0],
            [6, 1.0, 0.0, 1.0],
            [7, 1.0, 1.0, 1.0],
            [8, 0.0, 1.0, 1.0],
        ]
    )


def sample_model():
    model = FEModel(title="Test model", source="test")

    model.add_nodes(1, 0.0, 0.0, 0.0)
    model.add_nodes(2, 1.0, 0.0, 0.0)
    model.add_nodes(3, 0.0, 1.0, 0.0)
    model.add_nodes(4, 0.0, 0.0, 1.0)
    model.add_nodes(5, 1.0, 1.0, 1.0)

    model.add_elements(1, 1, [1, 2, 3, 4])
    model.add_elements(2, 2, [2, 3, 4, 5])
    model.add_part(part_id=1, name="part1", material_constants=[1, 1])
    model.add_part(part_id=2, name="part2", material_constants=[2, 1])

    model.section_info.append({"ID": 1, "constants": [1]})
    model.section_info.append({"ID": 2, "constants": [1]})

    model.material_info.append(
        {"type": "ELASTIC", "ID": 1, "constants": [1000.0, 210000.0, 0.3]}
    )

    model.material_info.append(
        {
            "type": "KELVIN-MAXWELL_VISCOELASTIC",
            "ID": 2,
            "constants": [1000.0, 210000.0, 0.3, 1.0, 1.0, 1.0, 1.0],
        }
    )

    return model


class TestFEModel:
    def test_femodel(self):
        mdl = sample_model()
        print(mdl)
        print(mdl.node_table)
        print(mdl.element_table)

        # check
        assert mdl.node_table.shape == (5, 4)
        assert mdl.element_table.shape == (2, 6)

        mdl.update_centroids()
        print(mdl.centroid_table)
        assert np.array(mdl.centroid_table).shape == (2, 3)

        assert len(mdl.part_info) == 2
        assert len(mdl.section_info) == 2
        assert len(mdl.material_info) == 2

    def test_add_nodes_single(self):
        # complete individual input
        mdl = sample_model()
        mdl.add_nodes(6, 2.0, 3.0, 4.0)

        np.testing.assert_equal(
            mdl.node_table[-1, :], np.array([6, 2.0, 3.0, 4.0])
        )

        # incomplete individual input
        with pytest.raises(ValueError):
            mdl.add_nodes(node_id=7, x=1.0, z=1.0)

        # repeat node ID
        with pytest.raises(ValueError):
            mdl.add_nodes(node_id=1, x=1.0, y=2.0, z=3.0)

        # complete array input
        mdl.add_nodes(node_array=np.array([7, 1.0, 2.0, 3.0]))

        # repeat node ID array
        with pytest.raises(ValueError):
            mdl.add_nodes(node_array=np.array([1, 1.0, 1.0, 1.0]))

    def test_add_nodes_multiple(self):
        mdl = sample_model()
        # array input
        node_arr = np.array([[6, 5.0, 5.0, 5.0], [7, 7.0, 7.0, 7.0]])

        mdl.add_nodes(node_array=node_arr)
        print(mdl.node_table)
        assert mdl.node_table.shape == (7, 4)

        node_arr = np.array([[6, 5.0, 5.0, 5.0], [8, 7.0, 7.0, 7.0]])
        # array input with repeat value
        with pytest.raises(ValueError):
            mdl.add_nodes(node_array=node_arr)

    def test_add_elements(self):
        # complete individual input
        mdl = sample_model()
        mdl.add_elements(3, 1, [1, 2, 3, 4])

        np.testing.assert_equal(
            mdl.element_table[-1, :],
            np.array([3, 1, 1.0, 2.0, 3.0, 4.0]),
        )

        # incomplete individual input
        with pytest.raises(ValueError):
            mdl.add_elements(element_id=3, nodes=[1, 2, 3, 4])

        # repeat node ID
        with pytest.raises(ValueError):
            mdl.add_elements(
                element_id=1,
                part_id=1,
                nodes=[
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
            )

        # complete array input
        mdl.add_elements(
            element_array=np.array([[5, 1, 1, 2, 3, 4], [6, 2, 1, 2, 3, 4]])
        )

        # repeat node ID array
        with pytest.raises(ValueError):
            mdl.add_elements(
                element_array=np.array(
                    [[1, 1, 1, 2, 3, 4], [7, 1, 1, 2, 3, 4]]
                )
            )

    def test_write_lsdyna_creates_file(self, tmp_path):
        mdl = sample_model()
        out_file = tmp_path / "test.k"
        mdl.write_lsdyna(os.path.join(tmp_path, "test.k"))
        assert out_file.exists()

        content = out_file.read_text()
        assert "*KEYWORD" in content
        assert "*NODE" in content
        assert "*ELEMENT_SOLID" in content
        assert "*PART" in content
        assert "*SECTION_SOLID" in content
        assert "*MAT_ELASTIC" in content
