from MRI2FE import FEModel

# Create model
model = FEModel(title="Test Model", source="test")

# Add sample data
model.add_node(1, 0.0, 0.0, 0.0)
model.add_node(2, 1.0, 0.0, 0.0)
model.add_node(3, 0.0, 1.0, 0.0)
model.add_node(4, 0.0, 0.0, 1.0)
model.add_node(5, 1.0, 1.0, 1.0)

model.add_element(1, [1, 2, 3, 4], part_id=1)
model.add_element(2, [2, 3, 4, 5], part_id=2)
model.add_part(part_id=1, name="part1", material_constants=[1, 1])
model.add_part(part_id=2, name="part2", material_constants=[2, 1])

model.section_info.append({"ID": 1, "constants": [1]})
model.section_info.append({"ID": 2, "constants": [1]})

model.material_info.append(
    {"type": "ELASTIC", "ID": 1, "constants": [1000.0, 210000.0, 0.3]}
)

model.material_info.append(
    {
        "type": "KELVIN-MAXWELL_VISCOELASTIC",
        "ID": 2,
        "constants": [1000.0, 210000.0, 0.3, 1.0, 1.0, 1.0, 1.0],
    }
)

# Write the LS-DYNA .k file
model.write_lsdyna("test/test_data/output_model.k")
