from MRI2FE.Pipelines.new_model import NewModel
import os

def test_mre_pipeline():
    model = NewModel()
    
    stiffness_path = "/Users/anshulshirude/Northeastern/Research/MRI2FE/src/MRI2FE/MRE/001Stiffness_warped.nii"
    damping_ratio_path = "/Users/anshulshirude/Northeastern/Research/MRI2FE/src/MRI2FE/MRE/001DR_warped.nii"
    mesh_path = "/Users/anshulshirude/Northeastern/Research/MRI2FE/src/MRI2FE/MRE/BrainWeb_Subject04_mesh.k"
    
    for file_path in [stiffness_path, damping_ratio_path, mesh_path]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
    
    print("Step 1: MRE/MRI Registration")
    registered = model.mre_mri_registration(
        mri_geometry_path=stiffness_path,
        mre_type="stiffness_damping",
        coregister_mre=True,
        stiffness_path=stiffness_path,
        damping_ratio_path=damping_ratio_path
    )
    print("Registration complete")
    
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
    
    print("\nStep 3: Calculate Element Centroids")
    mesh_data = model.calculate_element_centroids(
        mesh_path=mesh_path
    )
    print(f"Processed {len(mesh_data['centroids'])} elements")
    
    print("\nStep 4: Assign Material Properties")
    fe_model = model.assign_material_properties(
        material_constants=material_constants,
        mesh_data=mesh_data,
        mre_map=registered["mu"],
        output_csv="output_mapping.csv"
    )
    print("Material properties assigned to mesh")
    
    return fe_model

if __name__ == "__main__":
    try:
        fe_model = test_mre_pipeline()
        print("\nPipeline completed successfully!")
    except Exception as e:
        print(f"\nError in pipeline: {str(e)}")
