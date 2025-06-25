import femodel

# Create model
model = femodel.FEModel(title="Test Model", source="test")

# Add sample data
model.add_node(1, 0.0, 0.0, 0.0)
model.add_node(2, 1.0, 0.0, 0.0)
model.add_node(3, 0.0, 1.0, 0.0)
model.add_node(4, 0.0, 0.0, 1.0)

model.add_element(1, [1, 2, 3, 4], part_id=1)
model.add_part(1, {0: 1, 1: 1})

model.section_info.append({"ID": 1, "constants": [1]})

model.material_info.append({
    "type": "ELASTIC",
    "ID": 1,
    "constants": [210000, 0.3]
})

# Write the LS-DYNA .k file
model.write_lsdyna("output_model.k")
