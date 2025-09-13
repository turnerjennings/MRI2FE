from typing import List, Literal, Optional, Union, cast

import meshio
import numpy as np
from numpy.typing import ArrayLike


class FEModel:
    """Model object storing the generated FE model."""

    def __init__(
        self,
        title: str = "",
        source: str = "",
        imgout: Optional[str] = None,
        nodes: Optional[ArrayLike] = None,
        elements: Optional[ArrayLike] = None,
        parts: Optional[dict] = None,
        materials: Optional[List[dict]] = None,
        sections: Optional[List[dict]] = None,
    ):
        """Initialize the FE model

        Args:
            title (str, optional): Name for model, will be inserted as solver deck title on output. Defaults to "".
            source (str, optional): Source file for model, optional. Defaults to "".
            imgout (str, optional): Directory to save validation images to.  Defaults to None.
            nodes (Union[list, np.ndarray], optional): Array with shape (n,4) containing node ID, x, y, z coordinates. Defaults to None.
            elements (Union[list, np.ndarray], optional): Array with shape (n,m+2) containing element ID, part/group ID, and m connected nodes. Defaults to None.
            parts (dict, optional): Part definitions.  Dictionary keys are the part ID.  Each key is linked to a dictionary with entries "name" and "constants". Defaults to None.
            materials (List[dict], optional): Material definitions.  List of dictionaries.  Each dictionary must contain the keys "name", "ID", and "constants". Defaults to None.
            sections (List[dict], optional): Material section definitions.  List of dictionaries.  Each dictionary must contain the keys "ID" and "constants". Defaults to None.

        Raises:
            ValueError: Wrong input type.
        """
        self.metadata = {
            "title": title,
            "source": source,
            "num_nodes": 0,
            "num_elements": 0,
            "imgout": imgout,
        }

        self.centroid_table = None

        # create node table - array of shape (n,4) where each column has the format: [node_id, x, y, z]
        if nodes is not None:
            if isinstance(nodes, list):
                self.node_table: np.ndarray = np.array(nodes)
            elif isinstance(nodes, np.ndarray):
                self.node_table = nodes
            else:
                raise ValueError("Nodes must be a list or numpy array")

            self.metadata["num_nodes"] = self.node_table.shape[0]
        else:
            self.node_table = np.array([])
            self.metadata["num_nodes"] = 0

        # create element table - array of shape (n,m+2) where m is the number of nodes in the element type: [element_id, part_id, node1, node2, node3, node4...]
        if elements is not None:
            if isinstance(elements, np.ndarray):
                self.element_table = elements
            elif isinstance(elements, list):
                self.element_table = np.array(elements)
            else:
                raise ValueError("Elements must be a list or numpy array")
            self.metadata["num_elements"] = self.element_table.shape[0]
        else:
            self.element_table = np.array([])
            self.metadata["num_elements"] = 0

        # create centroid table - List of centroids: [x,y,z]
        self.centroid_table = None

        # create part info - Dictionary with keys "part id" and dictionary of "name":str and "constants":list
        if parts is not None:
            if isinstance(parts, dict):
                self.part_info = parts
            else:
                raise ValueError("Parts must be a dictionary")
        else:
            self.part_info = {}

        # create material info - List of Dictionaries with three entries: "type":str, "ID":int, and "constants":list[int,float]
        if materials is not None:
            if isinstance(materials, list):
                self.material_info = materials
            else:
                raise ValueError("Materials must be a list of dictionaries")
        else:
            self.material_info = []

        # create section info - List of Dictionaries with two entries: "ID": str and "constants":list[int, float]
        if sections is not None:
            if isinstance(sections, list):
                self.section_info = sections
            else:
                raise ValueError("Sections must be a list of dictionaries")
        else:
            self.section_info = []

    def add_nodes(
        self,
        node_id: Optional[int] = None,
        x: Optional[float] = None,
        y: Optional[float] = None,
        z: Optional[float] = None,
        node_array: Optional[ArrayLike] = None,
        force_insert: bool = False,
    ):
        """Add a node to the node table

        Args:
            node_id (int, optional): ID of the new node to add, must provide x-, y- and z- coordinates if specified. Defaults to None.
            x (float, optional): x-coordinate of the node, must also provide y- and z-coordinates. Defaults to None.
            y (float, optional): y-coordinate of the node, must also provide x- and z-coordinates. Defaults to None.
            z (float, optional): z-coordinate of the node, must also provide x- and y-coordinates. Defaults to None.
            node_array (np.ndarray, optional): (4,) array containing node_id,x,y,z coordinates. Defaults to None.
            force_insert (bool, optional): Whether to overwrite existing node with the same ID. Defaults to False.

        Raises:
            ValueError: Wrong input type.
            ValueError: Wrong array size.
            ValueError: Node ID already exists and force_insert false.
        """
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

        node_array = np.asarray(node_array)
        # check array dimensions match
        if (
            not indiv_input
            and node_array is not None
            and node_array.shape[1] != self.node_table.shape[1]
        ):
            raise ValueError(
                "Node array dimensions do not match node table dimensions"
            )

        if not force_insert and self.node_table.size > 0:
            # check if node already in the table
            node_table = np.atleast_2d(self.node_table)
            node_table = node_table[:, 0]
            matches = np.intersect1d(node_id_list, node_table)

            if len(matches) > 0:
                raise ValueError(
                    f"The following inserted nodes have duplicate IDs: {matches}"
                )

        # if node array is none, inserted nodes becomes array
        if self.node_table.size == 0 and indiv_input:
            self.node_table = np.array([node_id, x, y, z])
            self.metadata["num_nodes"] = 1

        elif self.node_table.size == 0:
            self.node_table = node_array
            self.metadata["num_nodes"] = node_array.shape[0]

        # add nodes to node table
        elif indiv_input:
            self.node_table = np.row_stack(
                (self.node_table, np.array([node_id, x, y, z]))
            )

            current_count = cast(int, self.metadata.get("num_nodes", 0))
            self.metadata["num_nodes"] = current_count + 1
        else:
            assert node_array is not None
            self.node_table = np.row_stack((self.node_table, node_array))

            current_count = cast(int, self.metadata.get("num_nodes", 0))
            self.metadata["num_nodes"] = current_count + node_array.shape[0]

    def add_elements(
        self,
        element_id: Optional[int] = None,
        part_id: Optional[int] = None,
        nodes: Optional[list] = None,
        element_array: Optional[np.ndarray] = None,
        force_insert: bool = False,
    ):
        """Add an element to the element table.

        Args:
            element_id (int, optional): Element ID, must also provide part_id and nodes if not None. Defaults to None.
            part_id (int, optional): Connected part ID, must also provide element_id and nodes if not None. Defaults to None.
            nodes (list, optional): Connected node IDs, must also provide element_id and part_id. Defaults to None.
            element_array (np.ndarray, optional): Array of shape (n+2,) containing element_id, part_id, n connected node ids. Defaults to None.
            force_insert (bool, optional): Whether to overwrite if element already exists with the same ID. Defaults to False.

        Raises:
            ValueError: Wrong input type.
            ValueError: Input shape mismatch with element table.
            ValueError: Element ID already exists and force_insert is False.
        """

        if all(var is not None for var in [element_id, part_id, nodes]):
            indiv_input = True
            element_id_list = [element_id]

        elif element_array is not None:
            indiv_input = False

            element_array = np.atleast_2d(element_array)
            element_id_list = element_array[:, 0].tolist()

        else:
            raise ValueError(
                "Must provide either (element_id, part_id, nodes) or element_array"
            )

        # check if array dimensions match
        if (
            not indiv_input
            and element_array is not None
            and not element_array.shape[1] == self.element_table.shape[1]
        ):
            raise ValueError(
                "Element array dimensions do not match self.element_table dimensions"
            )

        # check if element already exists
        if not force_insert and self.element_table.size > 0:
            # check if node already in the table
            element_table = np.atleast_2d(self.element_table)
            element_table = element_table[:, 0]
            matches = np.intersect1d(np.array(element_id_list), element_table)

            if len(matches) > 0:
                raise ValueError(
                    f"The following inserted nodes have duplicate IDs: {matches}"
                )

        # if element array is none, inserted element becomes array
        if self.element_table.size == 0 and indiv_input:
            assert nodes is not None
            self.element_table = np.array([element_id, part_id] + nodes)
            self.metadata["num_elements"] = 1

        elif self.element_table.size == 0:
            assert element_array is not None

            self.element_table = element_array
            self.metadata["num_elements"] = element_array.shape[0]

        elif indiv_input and nodes is not None:
            row_insert = np.array([element_id, part_id] + nodes)
            self.element_table = np.row_stack((self.element_table, row_insert))

            current_count = cast(int, self.metadata.get("num_elements", 0))
            self.metadata["num_elements"] = current_count + 1
        elif element_array is not None:
            self.element_table = np.row_stack(
                (self.element_table, element_array)
            )

            current_count = cast(int, self.metadata.get("num_elements", 0))
            self.metadata["num_elements"] = (
                current_count + element_array.shape[0]
            )

    def update_centroids(self):
        """Update the centroid table with all elements in the element table."""
        if self.element_table.size > 0:
            node_dict = {
                int(self.node_table[i, 0]): i
                for i in range(len(self.node_table))
            }
            n_elements = len(self.element_table)
            centroids = np.zeros((n_elements, 3))

            for i, element in enumerate(self.element_table):
                node_ids = np.unique(element[2:])
                node_indices = [
                    node_dict[int(node_id)] for node_id in node_ids
                ]
                centroids[i] = np.mean(
                    self.node_table[node_indices, 1:], axis=0
                )
            self.centroid_table = centroids

    def add_part(self, part_id: int, name: str, material_constants: list):
        """Add part information (e.g., material constants)."""
        self.part_info[part_id] = {}
        self.part_info[part_id]["name"] = name
        self.part_info[part_id]["constants"] = material_constants

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
            filename (str): Output file path for the .k file.
        """
        with open(filename, "w") as f:
            f.write("*KEYWORD\n")
            f.write(f"*TITLE\n{self.metadata['title']}\n")

            # Write nodes
            f.write("*NODE\n")
            for row in self.node_table:
                f.write(
                    f"{int(row[0]):8d}{row[1]:16.6f}{row[2]:16.6f}{row[3]:16.6f}\n"
                )

            # Write elements
            f.write("*ELEMENT_SOLID\n")
            for row in self.element_table:
                f.write(
                    f"{int(row[0]):>8d}{int(row[1]):>8d}\n"
                )  # eid and part id

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

                if "name" in mat:
                    f.write(f"*MAT_{mat_type}_TITLE\n")
                    f.write(f"{mat['name']}\n")
                else:
                    f.write(f"*MAT_{mat_type}\n")
                for item in mat_insert:
                    # convert from numpy types
                    if hasattr(
                        item, "item"
                    ):  # numpy scalars have .item() method
                        item = item.item()

                    # check and write type
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

    def from_meshio(
        self,
        mesh: Union[meshio.Mesh, str],
        element_type: Literal["tetra"] = "tetra",
        region_names: Optional[List[str]] = None,
    ):
        """Append mesh data from a meshio mesh object or file to the FEModel, offsetting IDs and
        updating metadata.

        Args:
            mesh (meshio.Mesh, str): Mesh object or filepath to mesh file.
            element_type ("tetra", optional): Element type to import. Only "tetra" is supported.
            region_names (List[str], optional): Optional part names for each region. Defaults to None.

        Raises:
            ValueError: If the element type of the new mesh does not match the existing FEModel.
            ValueError: If no elements of the specified type are found in the mesh.
        """
        if isinstance(mesh, str):
            mesh = meshio.read(mesh)

        # Check element type compatibility
        if self.element_table is not None and self.element_table.size > 0:
            # Assume all elements in self.element_table are of the same type
            num_nodes_per_elem = self.element_table.shape[1] - 2
            if element_type == "tetra" and num_nodes_per_elem != 4:
                raise ValueError(
                    "Element type mismatch: existing elements are not tetrahedrons."
                )

        mesh_nodes = mesh.points
        n_points = mesh_nodes.shape[0]
        # Offset node IDs
        max_node_id = (
            int(np.max(self.node_table[:, 0]))
            if self.node_table is not None and self.node_table.size > 0
            else 0
        )
        nids = np.arange(1, n_points + 1) + max_node_id
        new_nodes = np.column_stack((nids, mesh_nodes))

        # Find an element type in cells
        connectivity_index = None
        node_connectivity = None
        for idx, item in enumerate(mesh.cells):
            if item.type == element_type:
                # Map 0-based mesh indices to our new node IDs
                node_connectivity = nids[item.data]
                connectivity_index = idx
                break
        if node_connectivity is None:
            raise ValueError(
                f"Element type {element_type} not found in mesh cells"
            )

        n_elements = node_connectivity.shape[0]
        max_elem_id = (
            int(np.max(self.element_table[:, 0]))
            if self.element_table is not None and self.element_table.size > 0
            else 0
        )
        eids = np.arange(1, n_elements + 1) + max_elem_id

        # Offset part IDs
        if "medit:ref" not in mesh.cell_data:
            raise ValueError(
                "Mesh must contain 'medit:ref' cell data for part IDs"
            )

        pid = mesh.cell_data["medit:ref"][connectivity_index]
        max_part_id = (
            max([int(k) for k in self.part_info.keys()])
            if self.part_info
            else 0
        )
        pid_offset = max_part_id

        new_pid = pid + pid_offset
        new_elements = np.column_stack((eids, new_pid, node_connectivity))

        # Filter pid zero
        mask = new_pid > 0
        new_elements = new_elements[mask, :]
        new_pid_filtered = new_pid[mask]

        # all nodes are used
        new_nodes = new_nodes

        # Update part_info - only for parts that have elements (after filtering)
        for idx, id in enumerate(np.unique(new_pid_filtered)):
            if region_names is not None and idx < len(region_names):
                self.part_info[str(id)] = {
                    "name": region_names[idx],
                    "constants": [],
                }
            else:
                self.part_info[str(id)] = {"name": str(id), "constants": []}

        # Link node and element tables
        if self.node_table is not None and self.node_table.size > 0:
            self.node_table = np.row_stack((self.node_table, new_nodes))
        else:
            self.node_table = new_nodes

        if self.element_table is not None and self.element_table.size > 0:
            self.element_table = np.row_stack(
                (self.element_table, new_elements)
            )
        else:
            self.element_table = new_elements

        # Update metadata
        self.metadata["num_nodes"] = self.node_table.shape[0]
        self.metadata["num_elements"] = self.element_table.shape[0]
        if "source" in self.metadata and self.metadata["source"]:
            self.metadata["source"] = str(self.metadata["source"]) + (
                f", meshio:{getattr(mesh, 'filename', 'unknown')}"
            )
        else:
            self.metadata["source"] = (
                f"meshio:{getattr(mesh, 'filename', 'unknown')}"
            )
