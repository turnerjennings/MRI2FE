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
        cells={
            "hexahedron": [
                elem
                for part in model.elements.values()
                for elem in part["elements"]
            ]
        },
    )
    mesh.write(filename, file_format="abaqus")

    # Add more material types as needed underneath
    with open(filename, "a") as f:
        for mat_id, mat in model.materials.items():
            f.write(f"\n*MATERIAL, NAME={mat_id}\n")
            if mat["type"].lower() == "elastic":
                f.write("*ELASTIC\n")
                E = mat["properties"].get("E")
                nu = mat["properties"].get("nu")
                if E and nu:
                    f.write(f"{E}, {nu}\n")

        for part_id, part in model.elements.items():
            f.write(
                f"\n*SOLID SECTION, ELSET=part_{part_id}, MATERIAL={part['material']}\n"
            )
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

        material_id_map = {
            mat_id: idx
            for idx, mat_id in enumerate(model.materials.keys(), start=1)
        }

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
            MID = mat_idx

            if mat_type == "kelvin_maxwell":
                # Required parameters with error checking
                RO = props.get("density")
                E = props.get("E")
                PR = props.get("PR")
                if None in [RO, E, PR]:
                    raise ValueError(f"Missing required parameters for Kelvin-Maxwell material '{mat_id}'")

                # Parameters with defaults
                DC = props.get("DC", 0.0)
                FO = props.get("FO", 0.0)
                SO = props.get("SO", 0.0)

                # Compute derived parameters
                BULK = E / (3 * (1 - 2 * PR))
                G0 = props.get("G0", E / (2 * (1 + PR)))
                GI = props.get("GI", 0.0)

                f.write("\n*MAT_KELVIN_MAXWELL_VISCOELASTIC\n")
                f.write(f"{MID:8d}{RO:10.4e}{BULK:10.4e}{G0:10.4e}{GI:10.4e}{DC:10.4e}{FO:10.4e}{SO:10.4e}\n")

            elif mat_type == "elastic":
                # Required parameters with error checking
                RO = props.get("density")
                E = props.get("E")
                if RO is None or E is None:
                    raise ValueError(f"Missing RO or E for elastic material '{mat_id}'")

                # Parameters with defaults
                PR = props.get("PR", 0.0)
                DA = props.get("DA", 0.0)
                DB = props.get("DB", 0.0)
                K = props.get("K", 0.0)

                f.write("\n*MAT_ELASTIC\n")
                # Might need to handle negative E as integer
                f.write(
                    f"{MID:8d}{RO:10.4e}{E:10.4e}{PR:10.4e}{DA:10.4e}{DB:10.4e}{K:10.4e}\n"
                )

            else:
                raise NotImplementedError(f"Material type '{mat_type}' not supported.")

        f.write("\n*END\n")
