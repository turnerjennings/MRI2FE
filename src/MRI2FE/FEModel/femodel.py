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
        self.element_table = []  # List of elements: [element_id, part_id, node1, node2, node3, node4, ...]
        self.part_info = {}  # Dictionary of part_id -> material constants (e.g., Ginf, G1, Tau)

    def from_meshio(self, mesh: meshio.Mesh, default_part_id: int = 1):
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

        # Add nodes
        for node_id, (x, y, z) in enumerate(mesh.points, start=1):
            self.add_node(node_id, x, y, z)

        # Handle only one type of element for now (ex. "tetra")
        supported_keys = ["tetra"]
        found = False
        for key in supported_keys:
            if key in mesh.cells_dict:
                elements = mesh.cells_dict[key]
                for elem_id, node_ids in enumerate(elements, start=1):
                    # FIXED: + 1 offset since meshio is zero indexed but FEModel is one indexed
                    self.add_element(
                        elem_id, [i + 1 for i in node_ids], default_part_id
                    )
                found = True
                break

        if not found:
            raise ValueError(
                f"No supported cell types found in mesh. Supported: {supported_keys}"
            )

    def add_node(self, node_id: int, x: float, y: float, z: float):
        """Add a node to the node table."""
        self.node_table.append([node_id, x, y, z])
        self.metadata["num_nodes"] += 1

    def add_element(self, element_id: int, nodes: list, part_id: int):
        """Add an element to the element table.

        Args:
            element_id: The ID of the element
            nodes: List of node references
            part_id: The part ID for this element
        """
        self.element_table.append([element_id, part_id] + nodes)
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
