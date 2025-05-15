import pytest
from src.MRI2FE.output.fea_exporters import FEAModel, write_abaqus, write_lsdyna

def cube_nodes():
    return {
        1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)
    }


@pytest.fixture
def sample_model():
    return FEAModel(
        title="Test",
        nodes=cube_nodes(),
        elements={
            1: {
                "type": "C3D8",
                "elements": [[1, 2, 3, 4, 5, 6, 7, 8]],
                "material": "steel",
                "elform": 1,
                "aet": 0
            }
        },
        materials={
            "steel": {
                "type": "elastic",
                "properties": {"E": 210e9, "nu": 0.3, "density": 7850}
            }
        }
    )

@pytest.fixture
def visco_model():
    return FEAModel(
        title="Viscoelastic Test",
        nodes=cube_nodes(),
        elements={
            1: {
                "type": "C3D8",
                "elements": [[1, 2, 3, 4, 5, 6, 7, 8]],
                "material": "visco",
                "elform": 1,
                "aet": 0
            }
        },
        materials={
            "visco": {
                "type": "kelvin_maxwell",
                "properties": {
                    "E": 50e6,
                    "nu": 0.45,
                    "eta": 500.0,
                    "density": 1200
                }
            }
        }
    )

def test_write_abaqus_creates_file(tmp_path, sample_model):
    out_file = tmp_path / "test.inp"
    write_abaqus(sample_model, str(out_file))
    assert out_file.exists()

    content = out_file.read_text()
    assert "*NODE" in content or "*Node" in content
    assert "*ELEMENT" in content or "*Element" in content
    assert "*MATERIAL" in content
    assert "*SOLID SECTION" in content

def test_write_lsdyna_creates_file(tmp_path, sample_model):
    out_file = tmp_path / "test.k"
    write_lsdyna(sample_model, str(out_file))
    assert out_file.exists()

    content = out_file.read_text()
    assert "*KEYWORD" in content
    assert "*NODE" in content
    assert "*ELEMENT_SOLID" in content
    assert "*PART" in content
    assert "*SECTION_SOLID" in content
    assert "*MAT_ELASTIC" in content
    assert "7.8500e+03" in content or "7850.0000e+00" in content  # density value check
    assert "2.1000e+11" in content  # E value


def test_write_lsdyna_viscoelastic(tmp_path, visco_model):
    out_file = tmp_path / "visco.k"
    write_lsdyna(visco_model, str(out_file))
    assert out_file.exists()

    content = out_file.read_text()
    assert "*MAT_KELVIN_MAXWELL_VISCOELASTIC" in content
    assert "5.0000e+02" in content  # eta
    assert "1.2000e+03" in content or "1200.0000e+00" in content  # density
    assert "5.0000e+07" in content  # E
    assert "4.5000e-01" in content  # nu

