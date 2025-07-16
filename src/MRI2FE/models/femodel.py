import numpy as np
from typing import List, Union
from ..utilities import element_centroids
import meshio


class FEModel:
    def __init__(
        self,
        title: str = "",
        source: str = "",
        nodes: Union[list, np.ndarray] = None,
        elements: Union[list, np.ndarray] = None,
        parts: dict = None,
        materials: List[dict] = None,
        sections: List[dict] = None,
    ):
        """Initialize the FEModel data structure."""
        self.metadata = {
            "title": title,
            "source": source,
            "num_nodes": 0,
            "num_elements": 0,
        }

        # create node table - List of nodes: [node_id, x, y, z]
        if nodes is not None:
            if isinstance(nodes, np.ndarray):
                self.node_table = nodes.tolist()
            elif isinstance(nodes, np.ndarray):
                self.node_table = nodes
            else:
                raise ValueError("Nodes must be a list or numpy array")
        else:
            self.node_table: List[list] = []

        # create element table - List of elements: [element_id, part_id, node1, node2, node3, node4]
        if elements is not None:
            if isinstance(elements, np.ndarray):
                self.element_table = elements.tolist()
            elif isinstance(elements, list):
                self.element_table = elements
            else:
                raise ValueError("Elements must be a list or numpy array")

        else:
            self.element_table: List[list] = []

        # create centroid table - List of centroids: [x,y,z]
        self.centroid_table = []

        if len(self.element_table) > 0:
            for element in self.element_table:
                self.centroid_table.append(
                    element_centroids(element, np.array(self.node_table))
                )

        # create part info - Dictionary with keys "part id" and dictionary of "name":str and "constants":list
        if parts is not None:
            if isinstance(parts, dict):
                self.part_info = parts
            else:
                raise ValueError("Parts must be a dictionary")
        else:
            self.part_info: dict = {}

        # create material info - List of Dictionaries with three entries: "type":str, "ID":int, and "constants":list[int,float]
        if materials is not None:
            if isinstance(materials, list):
                self.material_info = materials
            else:
                raise ValueError("Materials must be a list of dictionaries")
        else:
            self.material_info: List[dict] = []

        # create section info - List of Dictionaries with two entries: "ID": str and "constants":list[int, float]
        if sections is not None:
            if isinstance(sections, list):
                self.section_info = sections
            else:
                raise ValueError("Sections must be a list of dictionaries")
        else:
            self.section_info: List[dict] = []

    def add_nodes(
        self,
        node_id: int = None,
        x: float = None,
        y: float = None,
        z: float = None,
        node_array: np.ndarray = None,
        force_insert: bool = False,
    ):
        """Add a node to the node table."""
        # check which input type is provided
        if all(var is not None for var in [node_id, x, y, z]):
            indiv_input = True
            node_id_list = np.array([node_id])

        elif node_array is not None:
            indiv_input = False
            node_array = np.atleast_2d(node_array)
            node_id_list = node_array[:, 0]

        else:
            raise ValueError(
                "Must provide either (node_id,x,y,z) or node_array"
            )

        if not force_insert and self.get_node_table().size > 0:
            # check if node already in the table
            node_table = self.get_node_table()
            node_table = np.atleast_2d(node_table)
            node_table = node_table[:, 0]
            matches = np.intersect1d(node_id_list, node_table)

            if len(matches) > 0:
                raise ValueError(
                    f"The following inserted nodes have duplicate IDs: {matches}"
                )

        # add nodes to node table
        if indiv_input:
            self.node_table.append([node_id, x, y, z])
            self.metadata["num_nodes"] += 1
        else:
            node_array = node_array.tolist()

            for node in node_array:
                self.node_table.append(node)
                self.metadata["num_nodes"] += 1

    def add_elements(
        self,
        element_id: int = None,
        part_id: int = None,
        nodes: list = None,
        element_array: np.ndarray = None,
        force_insert: bool = False,
    ):
        """Add an element to the element table.

        Args:
            element_id: The ID of the element
            nodes: List of node references
            part_id: The part ID for this element
        """

        if all(var is not None for var in [element_id, part_id, nodes]):
            indiv_input = True
            element_id_list = [element_id]

        elif element_array is not None:
            indiv_input = False

            element_array = np.atleast_2d(element_array)
            element_id_list = element_array[:, 0]

        else:
            raise ValueError(
                "Must provide either (element_id, part_id, nodes) or element_array"
            )

        # check if element already exists
        if not force_insert and self.get_element_table().size > 0:
            # check if node already in the table
            element_table = self.get_element_table()
            element_table = np.atleast_2d(element_table)
            element_table = element_table[:, 0]
            matches = np.intersect1d(element_id_list, element_table)

            if len(matches) > 0:
                raise ValueError(
                    f"The following inserted nodes have duplicate IDs: {matches}"
                )

        if indiv_input:
            self.element_table.append([element_id, part_id] + nodes)
            self.metadata["num_elements"] += 1
        else:
            element_list = element_array.tolist()

            for element in element_list:
                self.element_table.append(element)
                self.metadata["num_elements"] += 1

    def update_centroids(self):
        if self.get_element_table().size > 0:
            self.centroid_table = np.apply_along_axis(
                element_centroids,
                1,
                np.array(np.atleast_2d(self.element_table)),
                np.array(np.atleast_2d(self.node_table)),
            )
            self.centroid_table = self.centroid_table.tolist()

    def add_part(self, part_id: int, name: str, material_constants: list):
        """Add part information (e.g., material constants)."""
        self.part_info[part_id] = {}
        self.part_info[part_id]["name"] = name
        self.part_info[part_id]["constants"] = material_constants

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
            f.write(f"*TITLE\n{self.metadata['title']}\n")

            # Write nodes
            f.write("*NODE\n")
            for row in self.node_table:
                f.write(
                    f"{row[0]:8d}{row[1]:16.6f}{row[2]:16.6f}{row[3]:16.6f}\n"
                )

            # Write elements
            f.write("*ELEMENT_SOLID\n")
            for row in self.element_table:
                f.write(f"{row[0]:>8d}{row[1]:>8d}\n")  # eid and part id

                # write element connectivity, padding to 10-node format
                for i in range(2, len(row)):
                    f.write(f"{row[i]:>8d}")

                # Pad with last valid node or a dummy valid node ID (ex. repeat last node)
                last_node = row[-1]
                for i in range(len(row), 10):
                    f.write(f"{last_node:>8d}")

                # zero padding for n9 and n10
                f.write(f"{0:>8d}{0:>8d}")

                f.write("\n")

            # Writing parts
            for id, part in self.part_info.items():
                f.write("*PART\n")
                f.write(part["name"] + "\n")
                part_insert = [0, 0, 0, 0, 0, 0, 0, 0]
                part_insert[0] = int(id)
                # update default part with all available information
                for idx, item in enumerate(part["constants"]):
                    part_insert[idx + 1] = item

                for item in part_insert:
                    f.write(f"{item:>10d}")
                f.write("\n")

            # Write solid sections
            for sec in self.section_info:
                f.write("*SECTION_SOLID\n")
                secid = sec["ID"]
                elform = sec["constants"][0]
                aet = 0
                if len(sec["constants"]) > 1:
                    aet = sec["constants"][1]
                f.write(
                    f"{secid:>10d}{elform:>10d}{aet:10d}{0.0:>40.1f}{0.0:>10.1f}\n"
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

                mat_insert = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                mat_insert[0] = mat_id
                for idx, item in enumerate(props):
                    mat_insert[idx + 1] = item

                f.write(f"*MAT_{mat_type}\n")
                for item in mat_insert:
                    if isinstance(item, int):
                        f.write(f"{item:>10d}")
                    elif isinstance(item, float) and item == 0.0:
                        f.write(f"{0.0:>10.1f}")
                    elif isinstance(item, float):
                        f.write(f"{item:>10.2E}")
                    else:
                        raise ValueError(
                            f"Unsupported type in material constants: {type(item)}"
                        )

                f.write("\n")

            # End of file
            f.write("*END\n")


def model_from_meshio(
    mesh: Union[meshio.Mesh, str], title: str = "", source: str = ""
) -> FEModel:
    """
    Convert a meshio.Mesh object into a custom FEModel object.

    Args:
        mesh: meshio.Mesh object
        title: Metadata title for FEModel
        source: Metadata source for FEModel
        default_part_id: Default part ID for all elements

    Returns:
        FEModel instance with custom nodes and elements
    """
    if isinstance(mesh, str):
        mesh = meshio.read(mesh)

    # extract meshio points and create node numbering
    mesh_nodes = mesh.points
    n_points = mesh.points.shape[0]

    # apply offsets
    nids = np.arange(1, n_points + 1)
    new_nodes = np.column_stack((nids, mesh_nodes))

    # Handle only one type of element for now (ex. "tetra")
    supported_keys = ["tetra"]
    found = False
    for key in supported_keys:
        if key in mesh.cells_dict:
            # extract meshio cells and create node numbering

            mesh_elements = mesh.cells_dict[key]
            n_elems = mesh_elements.shape[0]

            # apply offsets
            eids = np.arange(1, n_elems + 1)
            mesh_elements = mesh_elements + 1

            new_elements = np.column_stack((eids, mesh_elements))

            found = True
            break

    if not found:
        raise ValueError(
            f"No supported cell types found in mesh. Supported: {supported_keys}"
        )

    out_mesh = FEModel(
        title=title, source=source, nodes=new_nodes, elements=new_elements
    )

    return out_mesh
