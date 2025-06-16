from MRI2FE.Pipelines.new_model import NewModel
from create_test_data import create_all_test_files
import os
import numpy as np
from ants.core.ants_image import ANTsImage

def validate_registration(registered):
    """Validate the registration step outputs"""
    assert "mu" in registered, "Missing stiffness (mu) in registered images"
    assert "xi" in registered, "Missing damping ratio (xi) in registered images"
    
    assert isinstance(registered["mu"], ANTsImage), "mu should be an ANTsImage"
    assert isinstance(registered["xi"], ANTsImage), "xi should be an ANTsImage"
    
    assert registered["mu"].shape == registered["xi"].shape, "mu and xi should have same dimensions"
    
    mu_array = registered["mu"].numpy()
    xi_array = registered["xi"].numpy()
    assert np.all(mu_array >= 0), "Stiffness values should be non-negative"
    assert np.all(xi_array >= 0), "Damping ratio values should be non-negative"
    
    print("✓ Registration validation passed")

def validate_material_constants(material_constants):
    """Validate the material constants calculation"""
    required_keys = ["Ginf", "G1", "Tau"]
    for key in required_keys:
        assert key in material_constants, f"Missing {key} in material constants"
    
    assert material_constants["Ginf"] > 0, "Ginf should be positive"
    assert material_constants["G1"] > 0, "G1 should be positive"
    assert material_constants["Tau"] > 0, "Tau should be positive"
    
    assert material_constants["Ginf"] < 1000, "Ginf seems unreasonably large"
    assert material_constants["G1"] < 1000, "G1 seems unreasonably large"
    assert material_constants["Tau"] < 100, "Tau seems unreasonably large"
    
    print("✓ Material constants validation passed")

def validate_mesh_data(mesh_data):
    """Validate the mesh data and centroids"""
    required_keys = ["element_table", "node_table", "centroids"]
    for key in required_keys:
        assert key in mesh_data, f"Missing {key} in mesh data"
    
    assert mesh_data["element_table"].shape[1] == 12, "Element table should have 12 columns"
    assert mesh_data["node_table"].shape[1] == 4, "Node table should have 4 columns"
    assert mesh_data["centroids"].shape[1] == 3, "Centroids should have 3 columns"
    
    assert len(mesh_data["centroids"]) == len(mesh_data["element_table"]), \
        "Number of centroids should match number of elements"
    
    min_coords = np.min(mesh_data["node_table"][:, 1:], axis=0)
    max_coords = np.max(mesh_data["node_table"][:, 1:], axis=0)
    assert np.all(mesh_data["centroids"] >= min_coords), "Centroids should be within node bounds"
    assert np.all(mesh_data["centroids"] <= max_coords), "Centroids should be within node bounds"
    
    print("✓ Mesh data validation passed")

def validate_fe_model(fe_model):
    """Validate the final FE model"""
    assert fe_model is not None, "FE model should not be None"
    
    node_table = fe_model.get_node_table()
    element_table = fe_model.get_element_table()
    
    assert len(node_table) > 0, "FE model should have nodes"
    assert len(element_table) > 0, "FE model should have elements"
    
    max_node_id = np.max(node_table[:, 0])
    print(f"\nDebug: Max node ID in node table: {max_node_id}")
    print(f"Debug: Number of nodes: {len(node_table)}")
    print(f"Debug: Number of elements: {len(element_table)}")
    
    element_table = element_table[element_table[:, 0].argsort()]
    
    for element in element_table:
        element_id = int(element[0])
        node_refs = element[2:]
        is_tet4 = all(node == 0 for node in node_refs[4:])
        active_nodes = node_refs[:4] if is_tet4 else node_refs
        
        assert np.all(active_nodes > 0), "Active node references should be positive"
        assert np.all(active_nodes <= max_node_id), "Active node references should be valid"
        
        if is_tet4:
            assert np.all(node_refs[4:] == 0), "Extra node references should be zero for TET4 elements"
    
    print("✓ FE model validation passed")

def test_mre_pipeline():
    stiffness_path, damping_ratio_path, mesh_path = create_all_test_files()
    model = NewModel()
    
    print("Step 1: MRE/MRI Registration")
    registered = model.mre_mri_registration(
        mri_geometry_path=stiffness_path,
        mre_type="stiffness_damping",
        coregister_mre=True,
        stiffness_path=stiffness_path,
        damping_ratio_path=damping_ratio_path
    )
    print("Registration complete")
    validate_registration(registered)
    
    print("\nStep 2: Calculate Material Constants")
    material_constants = model.calculate_material_constants(
        registered_images=registered,
        mre_type="stiffness_damping",
        mre_frequency=50.0
    )
    print("Material constants calculated:")
    print(f"Ginf: {material_constants['Ginf']}")
    print(f"G1: {material_constants['G1']}")
    print(f"Tau: {material_constants['Tau']}")
    validate_material_constants(material_constants)
    
    print("\nStep 3: Calculate Element Centroids")
    mesh_data = model.calculate_element_centroids(
        mesh_path=mesh_path
    )
    print(f"Processed {len(mesh_data['centroids'])} elements")
    validate_mesh_data(mesh_data)
    
    print("\nStep 4: Assign Material Properties")
    fe_model = model.assign_material_properties(
        material_constants=material_constants,
        mesh_data=mesh_data,
        mre_map=registered["mu"],
        output_csv="test_output_mapping.csv"
    )
    print("Material properties assigned to mesh")
    validate_fe_model(fe_model)
    
    return fe_model

if __name__ == "__main__":
    try:
        fe_model = test_mre_pipeline()
        print("\nPipeline completed successfully!")
    except AssertionError as e:
        print(f"\nValidation failed: {str(e)}")
    except Exception as e:
        print(f"\nError in pipeline: {str(e)}") 