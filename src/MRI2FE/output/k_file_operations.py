from ..FEModel.femodel import FEModel

import numpy as np
import os


def parse_k_file(fpath: str):
    """Load FE model from file and extract arrays for nodes and elements

    Args:
        fpath (str): filepath to the model to load

    Raises:
        TypeError: If fpath is not a string
        FileNotFoundError: If the file does not exist
        ValueError: If the file is empty or not in the expected format

    Returns:
        ect_array (np.array): (12,n) array of PID, EID, 10 node numbers
        node_array (np.array): (4,n) array of NID, xyz coordinates
    """
    # Validate input type
    if not isinstance(fpath, str):
        raise TypeError("fpath must be a string")

    # Check if file exists
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"File not found: {fpath}")

    # Check if file is empty
    if os.path.getsize(fpath) == 0:
        raise ValueError("File is empty")

    # load .k file
    try:
        with open(fpath, "r") as file:
            fcontent = file.read()
    except Exception as e:
        raise ValueError(f"Failed to read file: {str(e)}")

    lines = fcontent.split("\n")

    if not any(line.startswith("$#   eid     pid") for line in lines):
        raise ValueError("Invalid file format: Missing element data header")
    if "*NODE" not in fcontent:
        raise ValueError("Invalid file format: Missing node data section")

    element_section = fcontent.split("*ELEMENT_SOLID\n")[-1].split("*")[0]
    element_lines = [
        line.strip()
        for line in element_section.split("\n")
        if line.strip() and not line.startswith("$")
    ]

    data = []
    for line in element_lines:
        try:
            values = list(map(int, line.split()))
            if len(values) >= 10:
                if len(values) < 12:
                    values.extend([0] * (12 - len(values)))
                data.append(values)
        except ValueError:
            continue

    if not data:
        raise ValueError("No valid element data found in file")

    ect_array = np.array(data)

    node_section = fcontent.split("*NODE\n")[-1].split("*")[0]
    node_lines = [
        line.strip()
        for line in node_section.split("\n")
        if line.strip() and not line.startswith("$")
    ]

    node_data = []
    for line in node_lines:
        try:
            values = list(map(float, line.split()))
            if len(values) == 4:
                node_data.append(values)
        except ValueError:
            continue

    if not node_data:
        raise ValueError("No valid node data found in file")

    node_array = np.array(node_data)

    return ect_array, node_array


def write_head_k_file(fe_model: FEModel, fpath: str):
    """Write nodes and ect to an output file

    Args:
        ect (np.ndarray): Element connectivity table in LS-DYNA 10-node format
        nodecoords (np.ndarray): node ids and xyz coordinates
        fpath (str): Filepath to write output file to

    Raises:
        TypeError: If inputs are not of correct type
        ValueError: If inputs are None or if output directory doesn't exist
    """
    # Validate input types and None
    if fe_model is None:
        raise ValueError("fe_model cannot be None")
    if not isinstance(fe_model, FEModel):
        raise TypeError("fe_model must be a FEModel object")

    if fpath is None:
        raise ValueError("fpath cannot be None")
    if not isinstance(fpath, str):
        raise TypeError("fpath must be a string")

    # Validate output path exists
    output_dir = os.path.dirname(fpath)
    if output_dir and not os.path.exists(output_dir):
        raise ValueError(f"Output directory does not exist: {output_dir}")

    node_table = fe_model.get_node_table()
    element_table = fe_model.get_element_table()

    # format inputs
    k_file_boilerplate = [
        "*KEYWORD",
        "*ELEMENT_SOLID_TET4TOTET10",
        "$#   eid     pid",
        "$#    n1      n2      n3      n4      n5      n6      n7      n8      n9      n10",
        "*NODE",
        "#$ nid                x                y               z",
        "*END",
    ]

    eid_pid = element_table[:, 0:2]
    nodeconn = element_table[:, 2:]

    NID = node_table[:, 0]
    NXYZ = node_table[:, 1:]

    # define string templates
    eid_pid_widths = [8, 8]
    eid_pid_format_string = "".join(
        [f"{{:>{width}d}}" for width in eid_pid_widths]
    )

    nid_widths = [8, 8, 8, 8, 8, 8, 8, 8]
    nid_format_string = "".join([f"{{:>{width}d}}" for width in nid_widths])

    node_widths = [17, 17, 16]
    node_format_string = "".join([f"{{:>{width}f}}" for width in node_widths])

    # open the output file
    with open(fpath, "w") as file:
        # write starting boilerplate
        file.write(k_file_boilerplate[0] + "\n")

        # write ect
        file.write(k_file_boilerplate[1] + "\n")
        file.write(k_file_boilerplate[2] + "\n")
        file.write(eid_pid_format_string.format(*eid_pid[0, :]) + "\n")
        file.write(k_file_boilerplate[3] + "\n")
        file.write(nid_format_string.format(*nodeconn[0, :]) + "\n")

        for i in range(1, element_table.shape[0]):
            file.write(eid_pid_format_string.format(*eid_pid[i, :]) + "\n")
            file.write(nid_format_string.format(*nodeconn[i, :]) + "\n")

        # write nodes
        file.write(k_file_boilerplate[4] + "\n")
        # file.write(k_file_boilerplate[6] + "\n")

        for i in range(node_table.shape[0]):
            file.write(
                f"{int(NID[i]):>6d}"
                + node_format_string.format(*NXYZ[i, :])
                + "\n"
            )

        # write footer
        file.write(k_file_boilerplate[6] + "\n")


def edit_control_keyword(
    template: str,
    fpath: str,
    matprops: dict,
    includepath: str,
    title: str = None,
):
    """Create control keyword from template using a dictionary provided and save output

    Args:
        template (string): filepath to control keyword template .k file
        fpath (string): filepath to output the completed control keyword
        matprops (dict): dictionary of prony series constants for each desired model (G_inf, G1, Tau)
        includepath (string): filepath to the k file for the model to include
        title (string, optional): title for .k file if desired

    """

    entries = {}
    entries["include"] = includepath

    if title is not None:
        entries["title"] = title
    else:
        entries["title"] = ""

    n_segs = len(matprops["G1"])

    for i in range(1, n_segs + 1):
        entries[f"brain{i}G0"] = f"{matprops['Ginf'][i - 1]:>10.6f}"
        entries[f"brain{i}G1"] = f"{matprops['G1'][i - 1]:>10.6f}"
        entries[f"brain{i}ta"] = f"{matprops['Tau'][i - 1]:>10.6f}"

    for i in range(n_segs + 1, 10):
        entries[f"brain{i}G0"] = f"{0.0:>10.6f}"
        entries[f"brain{i}G1"] = f"{0.0:>10.6f}"
        entries[f"brain{i}ta"] = f"{0.0:>10.6f}"

    with open(template, "r") as fin:
        template = fin.read()

    template_filled = template.format(**entries)

    with open(fpath, "w") as fout:
        fout.write(template_filled)
