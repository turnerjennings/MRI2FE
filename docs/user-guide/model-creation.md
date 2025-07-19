# Model Creation

MRI2FE provides a comprehensive framework for creating and managing finite element models. The core of this functionality is the `FEModel` class, which handles nodes, elements, parts, materials, and sections.

## Overview

The model creation workflow involves:

1. **Creating a FEModel instance** with basic metadata
2. **Adding nodes and elements** to define the geometry
3. **Defining parts** to group elements
4. **Specifying materials** for material properties
5. **Adding sections** for element properties
6. **Exporting** to various formats (LS-DYNA, etc.)

## Core Classes

### `FEModel`

The main class for representing finite element models.

```python
from MRI2FE.models.femodel import FEModel

femodel = FEModel(
    title="Brain Model",
    source="MRI2FE generated",
    nodes=None,
    elements=None,
    parts=None,
    materials=None,
    sections=None
)
```

#### Constructor Parameters

- **`title`** (str, optional): Model title/name
- **`source`** (str, optional): Source information
- **`nodes`** (Union[list, np.ndarray], optional): Node data array
- **`elements`** (Union[list, np.ndarray], optional): Element data array
- **`parts`** (dict, optional): Part definitions
- **`materials`** (List[dict], optional): Material definitions
- **`sections`** (List[dict], optional): Section definitions

## Data Structures

### Node Table

Nodes are stored as a 2D array with format `[node_id, x, y, z]`:

```python
# Example node data
nodes = np.array([
    [1, 0.0, 0.0, 0.0],  # Node 1 at origin
    [2, 1.0, 0.0, 0.0],  # Node 2 at (1,0,0)
    [3, 0.0, 1.0, 0.0],  # Node 3 at (0,1,0)
    [4, 0.0, 0.0, 1.0]   # Node 4 at (0,0,1)
])
```

### Element Table

Elements are stored as a 2D array with format `[element_id, part_id, node1, node2, node3, node4, ...]`:

```python
# Example tetrahedral element data
elements = np.array([
    [1, 1, 1, 2, 3, 4],  # Element 1, Part 1, nodes 1,2,3,4
    [2, 1, 2, 3, 4, 5],  # Element 2, Part 1, nodes 2,3,4,5
])
```

### Part Information

Parts group elements and define material assignments:

```python
parts = {
    1: {"name": "brain_tissue", "constants": [1.0, 0.5, 0.1]},
    2: {"name": "skull", "constants": [2.0, 0.3, 0.05]}
}
```

### Material Information

Materials define constitutive properties:

```python
materials = [
    {"type": "viscoelastic", "ID": 1, "constants": [1.0, 0.5, 0.1]},
    {"type": "elastic", "ID": 2, "constants": [2.0, 0.3]}
]
```

### Section Information

Sections define element properties:

```python
sections = [
    {"ID": "solid", "constants": [1.0]},
    {"ID": "shell", "constants": [0.1, 0.1]}
]
```

## Core Methods

### Adding Nodes

```python
# Add individual node
femodel.add_nodes(
    node_id=1,
    x=0.0, y=0.0, z=0.0
)

# Add multiple nodes
node_array = np.array([
    [2, 1.0, 0.0, 0.0],
    [3, 0.0, 1.0, 0.0],
    [4, 0.0, 0.0, 1.0]
])
femodel.add_nodes(node_array=node_array)
```

### Adding Elements

```python
# Add individual element
femodel.add_elements(
    element_id=1,
    part_id=1,
    nodes=[1, 2, 3, 4]
)

# Add multiple elements
element_array = np.array([
    [2, 1, 2, 3, 4, 5],
    [3, 1, 3, 4, 5, 6]
])
femodel.add_elements(element_array=element_array)
```

### Managing Parts

```python
# Add a part
femodel.add_part(
    part_id=1,
    name="brain_tissue",
    material_constants=[1.0, 0.5, 0.1]
)

# Get part information
part_info = femodel.get_part_info()
print(part_info)
```

### Updating Centroids

```python
# Calculate element centroids
femodel.update_centroids()

# Access centroid data
centroids = femodel.centroid_table
```

## Creating Models from Meshes

### From meshio Mesh

```python
from MRI2FE.models.femodel import model_from_meshio
import meshio

# Load mesh
mesh = meshio.read("brain_mesh.vtk")

# Convert to FEModel
femodel = model_from_meshio(
    mesh=mesh,
    title="Brain Model",
    source="MRI2FE"
)
```

### From Existing Data

```python
import numpy as np
from MRI2FE.models.femodel import FEModel

# Create from arrays
nodes = np.array([
    [1, 0.0, 0.0, 0.0],
    [2, 1.0, 0.0, 0.0],
    [3, 0.0, 1.0, 0.0],
    [4, 0.0, 0.0, 1.0]
])

elements = np.array([
    [1, 1, 1, 2, 3, 4]
])

parts = {1: {"name": "tissue", "constants": [1.0, 0.5, 0.1]}}

femodel = FEModel(
    title="Simple Model",
    source="Manual creation",
    nodes=nodes,
    elements=elements,
    parts=parts
)
```

## Working with MRE Data

### Material Property Assignment

```python
from MRI2FE.Pipelines.new_model import NewModel

pipeline = NewModel()

# Calculate material constants from MRE data
material_constants = pipeline.calculate_material_constants(
    registered_images=registered_data,
    mre_type="stiffness_damping",
    mre_frequency=50.0
)

# Assign properties to mesh
femodel = pipeline.assign_material_properties(
    material_constants=material_constants,
    mesh_data=mesh_data,
    mre_map=mre_image,
    offset=3
)
```

### Element-Centroid Mapping

```python
# Calculate element centroids for MRE mapping
mesh_data = pipeline.calculate_element_centroids("brain_mesh.k")

# This provides centroid coordinates for mapping MRE properties
centroids = mesh_data["centroids"]
```

## Model Export

### LS-DYNA Format

```python
# Export to LS-DYNA .k file
femodel.write_lsdyna("output_model.k")
```

### Control Keywords

```python
from MRI2FE.output.k_file_operations import edit_control_keyword

# Create control keyword file
edit_control_keyword(
    template="template.k",
    fpath="control.k",
    matprops=material_properties,
    includepath="model.k",
    title="Brain Impact Simulation"
)
```

## Examples

### Basic Model Creation

```python
import numpy as np
from MRI2FE.models.femodel import FEModel

def create_simple_tetrahedron():
    """Create a simple tetrahedral model."""

    # Define nodes
    nodes = np.array([
        [1, 0.0, 0.0, 0.0],
        [2, 1.0, 0.0, 0.0],
        [3, 0.0, 1.0, 0.0],
        [4, 0.0, 0.0, 1.0]
    ])

    # Define elements
    elements = np.array([
        [1, 1, 1, 2, 3, 4]
    ])

    # Define parts
    parts = {1: {"name": "tissue", "constants": [1.0, 0.5, 0.1]}}

    # Create model
    femodel = FEModel(
        title="Simple Tetrahedron",
        source="Example",
        nodes=nodes,
        elements=elements,
        parts=parts
    )

    return femodel

# Usage
model = create_simple_tetrahedron()
model.write_lsdyna("tetrahedron.k")
```

### Multi-Part Model

```python
def create_multi_part_model():
    """Create a model with multiple parts."""

    # Create model
    femodel = FEModel(title="Multi-Part Model")

    # Add nodes
    nodes = np.array([
        [1, 0.0, 0.0, 0.0],
        [2, 1.0, 0.0, 0.0],
        [3, 0.0, 1.0, 0.0],
        [4, 0.0, 0.0, 1.0],
        [5, 2.0, 0.0, 0.0],
        [6, 1.0, 1.0, 0.0],
        [7, 1.0, 0.0, 1.0],
        [8, 2.0, 1.0, 1.0]
    ])

    femodel.add_nodes(node_array=nodes)

    # Add elements for part 1 (brain)
    brain_elements = np.array([
        [1, 1, 1, 2, 3, 4],
        [2, 1, 2, 3, 4, 5]
    ])
    femodel.add_elements(element_array=brain_elements)

    # Add elements for part 2 (skull)
    skull_elements = np.array([
        [3, 2, 5, 6, 7, 8],
        [4, 2, 6, 7, 8, 1]
    ])
    femodel.add_elements(element_array=skull_elements)

    # Add parts
    femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])
    femodel.add_part(2, "skull", [2.0, 0.3, 0.05])

    return femodel
```

### Model from Mesh

```python
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio

def create_model_from_mri(mri_path, output_path):
    """Create FE model from MRI data."""

    # Generate mesh
    mesh = mesh_from_nifti(mri_path, optimize=True)

    # Convert to FE model
    femodel = model_from_meshio(
        mesh=mesh,
        title="Brain Model",
        source="MRI2FE"
    )

    # Add default part
    femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])

    # Save model
    femodel.write_lsdyna(output_path)

    return femodel
```

## Advanced Features

### Model Validation

```python
def validate_model(femodel):
    """Validate FE model for common issues."""

    issues = []

    # Check for duplicate nodes
    node_ids = femodel.node_table[:, 0]
    if len(node_ids) != len(set(node_ids)):
        issues.append("Duplicate node IDs found")

    # Check for duplicate elements
    element_ids = femodel.element_table[:, 0]
    if len(element_ids) != len(set(element_ids)):
        issues.append("Duplicate element IDs found")

    # Check element connectivity
    for i, element in enumerate(femodel.element_table):
        node_refs = element[2:]  # Skip element_id and part_id
        for node_ref in node_refs:
            if node_ref not in node_ids:
                issues.append(f"Element {i+1} references non-existent node {node_ref}")

    return issues

# Usage
model = create_simple_tetrahedron()
issues = validate_model(model)
if issues:
    print("Model validation issues:", issues)
else:
    print("Model is valid!")
```

### Model Statistics

```python
def print_model_stats(femodel):
    """Print comprehensive model statistics."""

    print(f"Model: {femodel.metadata['title']}")
    print(f"Source: {femodel.metadata['source']}")
    print(f"Nodes: {femodel.metadata['num_nodes']}")
    print(f"Elements: {femodel.metadata['num_elements']}")

    # Part statistics
    for part_id, part_info in femodel.part_info.items():
        elements_in_part = len(femodel.element_table[femodel.element_table[:, 1] == part_id])
        print(f"Part {part_id} ({part_info['name']}): {elements_in_part} elements")

    # Bounds
    if femodel.node_table is not None:
        coords = femodel.node_table[:, 1:]
        bounds = np.column_stack([np.min(coords, axis=0), np.max(coords, axis=0)])
        print(f"Bounds: X[{bounds[0,0]:.2f}, {bounds[0,1]:.2f}], "
              f"Y[{bounds[1,0]:.2f}, {bounds[1,1]:.2f}], "
              f"Z[{bounds[2,0]:.2f}, {bounds[2,1]:.2f}]")
```

## Troubleshooting

### Common Issues

1. **Node ID conflicts:**

   ```python
   # Use force_insert=True to overwrite existing nodes
   femodel.add_nodes(node_id=1, x=0, y=0, z=0, force_insert=True)
   ```

2. **Element connectivity errors:**

   ```python
   # Ensure all referenced nodes exist
   existing_node_ids = set(femodel.node_table[:, 0])
   for element in elements:
       for node_ref in element[2:]:
           if node_ref not in existing_node_ids:
               print(f"Warning: Node {node_ref} not found")
   ```

3. **Memory issues with large models:**
   ```python
   # Process in batches for large datasets
   batch_size = 1000
   for i in range(0, len(nodes), batch_size):
       batch = nodes[i:i+batch_size]
       femodel.add_nodes(node_array=batch)
   ```

## Best Practices

1. **Use consistent node/element numbering** (1-indexed recommended)
2. **Validate models** before export
3. **Use descriptive part names** for clarity
4. **Include proper metadata** (title, source)
5. **Test with small models** before scaling up
6. **Backup models** before major modifications

## References

- [LS-DYNA Keyword Manual](https://www.dynasupport.com/manuals)
- [Finite Element Method Fundamentals](https://en.wikipedia.org/wiki/Finite_element_method)
