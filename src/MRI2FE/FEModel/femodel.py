import numpy as np
import meshio
from typing import List, Tuple


class FEModel:
    def __init__(self, title: str = "", source: str = ""):
        """Initialize the FEModel data structure."""
        self.metadata = {
            "title": title,
            "source": source,
            "num_nodes": 0,
            "num_elements": 0,
        }
        self.node_table: List[
            Tuple[int, float, float, float]
        ] = []  # List of nodes: [node_id, x, y, z]

        self.element_table: List[
            Tuple[int, int, int, int, int, int]
        ] = []  # List of elements: [element_id, part_id, node1, node2, node3, node4]

        self.part_info: dict = {}  # Dictionary with keys "part id" and list of assigned constants [section id, material id]

        self.material_info: List[
            dict
        ] = []  # List of Dictionaries with three entries: "type":str, "ID":int, and "constants":list[int,float]
        
        self.section_info: List[
            dict
        ] = []  # List of Dictionaries with two entries: "ID": str and "constants":list[int, float]

    def add_node(self, node_id: int, x: float, y: float, z: float):
        """Add a node to the node table."""
        self.node_table.append([node_id, x, y, z])
        self.metadata["num_nodes"] += 1

    def add_element(self, element_id: int, nodes: list, part_id: int):
        """Add an element to the element table."""
        self.element_table.append([element_id] + nodes + [part_id])
        self.metadata["num_elements"] += 1

    def add_part(self, part_id: int, material_constants: dict):
        """Add part information (e.g., material constants)."""
        self.part_info[part_id] = material_constants

    def get_node_table(self):
        """Return the node table."""
        return np.array(self.node_table)

    def get_element_table(self):
        """Return the element table."""
        return np.array(self.element_table)

    def get_part_info(self):
        """Return the part information."""
        return self.part_info

    def __repr__(self):
        """String representation of the FEModel."""
        return (
            f"FEModel(title={self.metadata['title']}, "
            f"source={self.metadata['source']}, "
            f"num_nodes={self.metadata['num_nodes']}, "
            f"num_elements={self.metadata['num_elements']})"
        )

    def write_lsdyna(self, filename: str):
        """
        Write FE model data to an LS-DYNA .k file.

        Args:
            model (FEAModel): Model containing nodes, elements, and materials.
            filename (str): Output file path for the .k file.
        """
        with open(filename, "w") as f:
            f.write("*KEYWORD\n")
            f.write(f"*TITLE\n{self.metadata['title']}\n\n")

            # Write nodes
            f.write("*NODE\n")
            for row in self.node_table:
                f.write(
                    f"{row[0]:8d}{row[1]:16.6f}{row[2]:16.6f}{row[3]:16.6f}\n"
                )

            # Write elements
            for row in self.element_table:
                f.write("\n*ELEMENT_SOLID\n")
                f.write(f"{row[0]:>8d}{row[1]:>8d}\n")  # eid and part id

                # write element connectivity, padding to 10-node format
                for i in range(2, len(row)):
                    f.write(f"{row[i]:>8d}")

                for i in range(len(row), 11):
                    f.write(f"{0:>8d}")

                f.write("\n")

            # Writing parts
            f.write("\n*PART\n")
            for id, part in self.part_info.items():
                part_insert = [0, 0, 0, 0, 0, 0, 0, 0]
                part_insert[0] = int(id)
                # update default part with all available information
                for idx, item in enumerate(part):
                    part_insert[idx + 1] = item

                for item in part_insert:
                    f.write(f"{item:>8d}")
                f.write("\n")

            # Write solid sections
            for sec in self.section_info:
                f.write("\n*SECTION_SOLID\n")
                secid = sec["ID"]
                elform = sec["constants"][0]
                aet = 0
                if len(sec["constants"] > 1):
                    aet = sec["constants"][1]
                f.write(
                    f"{secid:>8d}{elform:>8d}{aet:8d}{0.0:>32.1f}{0.0:>8.1f}\n"
                )

            # Write materials
            for mat in self.material_info:
                mat_type = mat["type"]
                mat_id = mat["ID"]
                props = mat["constants"]

                # check if multi-line input card, raise error if yes
                if len(props) > 7:
                    raise ValueError(
                        f"Error in material id {mat_id}: multi-line input cards not supported"
                    )

                mat_insert = [0, 0, 0, 0, 0, 0, 0, 0]
                mat_insert[0] = mat_id
                for idx, item in enumerate(props):
                    mat_insert[idx + 1] = item

                f.write(f"*MAT_{mat_type}\n")
                for item in mat_insert:
                    f.write(f"{item:>8d}")

                f.write("\n")

            # End of file
            f.write("\n*END\n")

    def write_abaqus(self, filename: str):
        """
        Write FE model data to an ABAQUS .inp file.

        Args:
            model (FEAModel): Model containing nodes, elements, and materials.
            filename (str): Output file path for the .inp file.
        """

        nodes = self.get_node_table()

        cells = self.get_element_table()
        mesh = meshio.Mesh(
            points=nodes,
            cells={"tetra": cells},
        )
        mesh.write(filename, file_format="abaqus")

        # Add more material types as needed underneath
        with open(filename, "a") as f:
            for mat_id, mat in self.material_info.items():
                f.write(f"\n*MATERIAL, NAME={mat_id}\n")
                if mat["type"].lower() == "elastic":
                    f.write("*ELASTIC\n")
                    youngs_modulus = mat["properties"].get("E")
                    poisson_ratio = mat["properties"].get("nu")
                    if youngs_modulus and poisson_ratio:
                        f.write(f"{youngs_modulus}, {poisson_ratio}\n")

            for part_id, part in self.part_info.items():
                f.write(
                    f"\n*SOLID SECTION, ELSET=part_{part_id}, MATERIAL={part['material']}\n"
                )
                f.write("1.0\n")
