import pytest
from MRI2FE import FEModel
import numpy as np
import os

def cube_nodes() -> np.ndarray:
    return np.array([
        [1, 0., 0., 0.],
        [2, 1., 0., 0.],
        [3, 1., 1., 0.],
        [4, 0., 1., 0.],
        [5, 0., 0., 1.],
        [6, 1., 0., 1.],
        [7, 1., 1., 1.],
        [8, 0., 1., 1.],
    ])


@pytest.fixture
def sample_model():

    model = FEModel(title="Test model",
                    source="test")

    model.add_nodes(1, 0.0, 0.0, 0.0)
    model.add_nodes(2, 1.0, 0.0, 0.0)
    model.add_nodes(3, 0.0, 1.0, 0.0)
    model.add_nodes(4, 0.0, 0.0, 1.0)
    model.add_nodes(5, 1.0, 1.0, 1.0)

    model.add_elements(1,1, [1, 2, 3, 4])
    model.add_elements(2,2, [2, 3, 4, 5])
    model.add_part(part_id=1, name="part1", material_constants=[1,1])
    model.add_part(part_id=2, name="part2", material_constants=[2,1])

    model.section_info.append({"ID": 1, "constants": [1]})
    model.section_info.append({"ID": 2, "constants": [1]})

    model.material_info.append({
        "type": "ELASTIC",
        "ID": 1,
        "constants": [1000., 210000., 0.3]
    })

    model.material_info.append({
        "type": "KELVIN-MAXWELL_VISCOELASTIC",
        "ID": 2,
        "constants": [1000., 210000., 0.3, 1.0, 1.0, 1.0, 1.0]
    })

    return model


def test_write_lsdyna_creates_file(tmp_path, sample_model):
    out_file = tmp_path / "test.k"
    sample_model.write_lsdyna(os.path.join(tmp_path,"test.k"))
    assert out_file.exists()

    content = out_file.read_text()
    assert "*KEYWORD" in content
    assert "*NODE" in content
    assert "*ELEMENT_SOLID" in content
    assert "*PART" in content
    assert "*SECTION_SOLID" in content
    assert "*MAT_ELASTIC" in content
