import numpy as np


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
