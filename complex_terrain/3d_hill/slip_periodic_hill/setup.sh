# convenience script for setting up the ic/mesh
aprepro -c# periodic_hill.yaml run.yaml
abl_mesh -i run.yaml
python create_complex_terrain.py
nalu_preprocess -i run.yaml
