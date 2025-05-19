import numpy as np
import meshio

class FEModel:
    def __init__(self, title: str = "", source: str = ""):
        """Initialize the FEModel data structure."""
        self.metadata = {
            "title": title,
            "source": source,
            "num_nodes": 0,
            "num_elements": 0,
        }
        self.node_table = []  # List of nodes: [node_id, x, y, z]
        self.element_table = []  # List of elements: [element_id, node1, node2, node3, node4, part_id]
        self.part_info = {}  # Dictionary of part_id -> material constants (e.g., Ginf, G1, Tau)
        self.materials = [] # List of dictionaries with material assignments

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
    
    def write_abaqus(self, filename: str):
        """
        Write FE model data to an ABAQUS .inp file.
        
        Args:
            model (FEAModel): Model containing nodes, elements, and materials.
            filename (str): Output file path for the .inp file.
        """
        points=np.array([node[1:] for node in self.node_table]),

        tetra_cells = []
        part_data = []

        for elem in self.element_table:
            # get element node ids
            node_indices = elem[1:5]
            
            tetra_cells.append(node_indices)
            
            # get part ID
            part_id = elem[5]
            part_data.append(part_id)

        mesh = meshio.Mesh(
            points = points,
            cells = {"tetra": np.array(tetra_cells)},
            cell_data = {"part_id": [np.array(part_data)]}
        )
        mesh.write(filename, file_format="abaqus")

        # Add more material types as needed underneath
        with open(filename, "a") as f:
            for mat_id, mat in enumerate(self.materials):
                f.write(f"\n*MATERIAL, NAME={mat_id}\n")
                if mat["type"].lower() == "elastic":
                    f.write("*ELASTIC\n")
                    if "E" in mat and "nu" in mat:
                        youngs_modulus = mat["E"]
                        poisson_ratio = mat["nu"]

                        f.write(f"{youngs_modulus}, {poisson_ratio}\n")

        #write sections
        #TODO fix part/section ID data structure
        for part_id, part in self.elements.items():
            f.write(
                f"\n*SOLID SECTION, ELSET=part_{part_id}, MATERIAL={part['material']}\n"
            )
            f.write("1.0\n")

    def write_lsdyna(self, filename: str):
        """
        Write FE model data to an LS-DYNA .k file.

        Args:
            model (FEAModel): Model containing nodes, elements, and materials.
            filename (str): Output file path for the .k file.
        """
        with open(filename, "w") as f:
            f.write("*KEYWORD\n")
            if "title" in self.metadata:
                f.write(f"*TITLE\n{self.metadata["title"]}\n\n")
            else:
                f.write(f"*TITLE\nSolver Deck\n\n")

            # Write nodes
            f.write("*NODE\n")
            for n in self.node_table:

                id = n[0]
                x = n[1]
                y = n[2]
                z = n[3]

                f.write(f"{id:8d}{x:16.6f}{y:16.6f}{z:16.6f}\n")

            # Write elements
            f.write("\n*ELEMENT_SOLID\n")
            for element in self.element_table:

                eid = element[0]
                pid = element[5]

                nodes = element[1:5]

                while len(nodes) < 10:
                    nodes.append(0)
                
                f.write(f"{eid:8d}{pid:8d}\n")

                nodes = "".join(f"{n:8d}" for n in nodes)
                f.write(f"{nodes}\n")              


            material_id_map = {
                mat_id: idx
                for idx, mat_id in enumerate(model.materials.keys(), start=1)
            }

            #TODO fix data structure
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
            for mat_idx, (mat_id, mat) in enumerate(
                model.materials.items(), start=1
            ):
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
                    f.write(
                        f"{mat_idx:8d}{density:10.4e}{youngs_modulus:10.4e}{poisson_ratio:10.4e}{gk:10.4e}{viscosity_coefficient:10.4e}\n"
                    )

                elif mat_type == "elastic":
                    if youngs_modulus is None or poisson_ratio is None:
                        raise ValueError(f"Missing E or nu for material {mat_id}")
                    f.write("\n*MAT_ELASTIC\n")
                    f.write(
                        f"{mat_idx:8d}{density:10.4e}{youngs_modulus:10.4e}{poisson_ratio:10.4e}\n"
                    )

                else:
                    raise NotImplementedError(
                        f"Material type '{mat_type}' not supported."
                    )

            # End of file
            f.write("\n*END\n")
