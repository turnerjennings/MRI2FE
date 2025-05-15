import meshio
from dataclasses import dataclass

@dataclass
class FEAModel:
    """
    Container for finite element model data.

    Attributes:
        title (str): Model title or description.
        nodes (dict): {node_id: (x, y, z)} node coordinates.
        elements (dict): {part_id: {"type": str, "elements": list, "material": str}}
        materials (dict): {material_id: {"type": str, "properties": dict}}
    """
    title: str
    nodes: dict
    elements: dict
    materials: dict

def write_abaqus(model: FEAModel, filename: str):
    """
    Write FE model data to an ABAQUS .inp file.

    Args:
        model (FEAModel): Model containing nodes, elements, and materials.
        filename (str): Output file path for the .inp file.
    """
    mesh = meshio.Mesh(
        points=[coord for coord in model.nodes.values()],
        cells={"hexahedron": [
            elem for part in model.elements.values()
            for elem in part["elements"]
        ]}
    )
    mesh.write(filename, file_format="abaqus")

    # Add more material types as needed underneath
    with open(filename, "a") as f:
        for mat_id, mat in model.materials.items():
            f.write(f"\n*MATERIAL, NAME={mat_id}\n")
            if mat["type"].lower() == "elastic":
                f.write("*ELASTIC\n")
                youngs_modulus = mat["properties"].get("E")
                poisson_ratio = mat["properties"].get("nu")
                if youngs_modulus and poisson_ratio:
                    f.write(f"{youngs_modulus}, {poisson_ratio}\n")

        for part_id, part in model.elements.items():
            f.write(f"\n*SOLID SECTION, ELSET=part_{part_id}, MATERIAL={part['material']}\n")
            f.write("1.0\n")

def write_lsdyna(model: FEAModel, filename: str):
    """
    Write FE model data to an LS-DYNA .k file.

    Args:
        model (FEAModel): Model containing nodes, elements, and materials.
        filename (str): Output file path for the .k file.
    """
    with open(filename, "w") as f:
        f.write("*KEYWORD\n")
        f.write(f"*TITLE\n{model.title}\n\n")

        # Write nodes
        f.write("*NODE\n")
        for nid, (x, y, z) in model.nodes.items():
            f.write(f"{nid:8d}{x:16.6f}{y:16.6f}{z:16.6f}\n")

        # Write elements
        for part_id, part in model.elements.items():
            f.write("\n*ELEMENT_SOLID\n")
            for eid, elem in enumerate(part["elements"], start=1):
                nodes = "".join(f"{nid:8d}" for nid in elem)
                f.write(f"{eid:8d}{part_id:8d}{nodes}\n")

        material_id_map = {mat_id: idx for idx, mat_id in enumerate(model.materials.keys(), start=1)}

        # When writing parts
        f.write("\n*PART\n")
        for part_id, part in model.elements.items():
            mid = material_id_map[part["material"]]
            f.write(f"{part_id:8d}{part_id:8d}{mid:8d}\n")

        # Write solid sections
        f.write("\n*SECTION_SOLID\n")
        for part_id, part in model.elements.items():
            elform = part.get("elform", 1)
            aet = part.get("aet", 0)
            f.write(f"{part_id:8d}{elform:8d}{aet:8d}\n")

        # Write materials
        for mat_idx, (mat_id, mat) in enumerate(model.materials.items(), start=1):
            mat_type = mat["type"].lower()
            props = mat["properties"]

            youngs_modulus = props.get("E")
            poisson_ratio = props.get("nu")
            density = props.get("density", 1.0)

            if mat_type == "kelvin_maxwell":
                viscosity_coefficient = props.get("eta", 0.0)
                if youngs_modulus is None or poisson_ratio is None:
                    raise ValueError(f"Missing E or nu for material {mat_id}")
                gk = youngs_modulus / (3 * (1 - 2 * poisson_ratio))

                f.write("\n*MAT_KELVIN_MAXWELL_VISCOELASTIC\n")
                f.write(f"{mat_idx:8d}{density:10.4e}{youngs_modulus:10.4e}{poisson_ratio:10.4e}{gk:10.4e}{viscosity_coefficient:10.4e}\n")

            elif mat_type == "elastic":
                if youngs_modulus is None or poisson_ratio is None:
                    raise ValueError(f"Missing E or nu for material {mat_id}")
                f.write("\n*MAT_ELASTIC\n")
                f.write(f"{mat_idx:8d}{density:10.4e}{youngs_modulus:10.4e}{poisson_ratio:10.4e}\n")

            else:
                raise NotImplementedError(f"Material type '{mat_type}' not supported.")

        # End of file
        f.write("\n*END\n")