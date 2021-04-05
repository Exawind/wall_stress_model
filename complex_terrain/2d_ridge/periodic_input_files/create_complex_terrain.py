# -*- coding: utf-8 -*-

"""\
Basic usage of pySTK interface
==============================

This tutorial provides a high-level overview of using pySTK API to interact with
Exodus-II databases and the Sierra Toolkit library.

The main entry point is the StkMesh instance which is a Python-wrapper that
holds the MetaData, BulkData, and StkMeshIoBroker instances. The user can
create multiple instances of StkMesh to manipulate different meshes.

Every class in STK (e.g., MetaData, BulkData, Part, Selector, etc.) has a
python wrapped counterpart with the `Stk` prefix. For example, `MetaData` is
exposed as `StkMetaData` within the python layer

"""

import numpy as np
import stk
from stk import StkMesh, Parallel, StkState, StkSelector, StkRank

par = Parallel.initialize()
print("MPI: rank = %d; size = %d"%(par.rank, par.size))

# Create a STK mesh instance that holds MetaData and BulkData
mesh = StkMesh(par, ndim=3)


mesh.read_mesh_meta_data(
    "mesh.exo",
    purpose=stk.DatabasePurpose.READ_MESH,
    auto_decomp=True,
    auto_declare_fields=True,
    auto_decomp_type="rcb")

# Equivalent statement as above
# mesh.read_mesh_meta_data("periodic3d.g")

# At this point the mesh metadata has been read in, but the bulk data hasn't
# been populated yet. The metadata is open for modification, i.e., the user can
# add parts, fields, etc.

# Register fields
#disp = mesh.meta.declare_vector_field("mesh_displacement")

# Fields of type other than `double` can be declared with special syntax
iblank = mesh.meta.declare_scalar_field_t[int]("iblank")

# Register the fields on desired parts
#disp.add_to_part(mesh.meta.universal_part,
#                     mesh.meta.spatial_dimension)
iblank.add_to_part(mesh.meta.universal_part)

# Commit the metadata and load the mesh. Also create edges at this point (default is False)
mesh.populate_bulk_data(create_edges=True)

mesh.stkio.read_defined_input_fields(0.0)

print("Metadata is committed: ", mesh.meta.is_committed)

# Access the coordinate field
coords = mesh.meta.coordinate_field
print("Coordinates field: ", coords.name)

# Access a part by name
part = mesh.meta.get_part("fluid_part")
print("Part exists = ", (not part.is_null))

# Get a stk::mesh::Selector instance for all locally-owned entities of this part
sel = part & mesh.meta.locally_owned_part

# Check if the selector is empty for a particular entity type
print("Fluid_part has elems: ", not sel.is_empty(StkRank.ELEM_RANK))

###
### Common looping concepts
###

# Iterating over entities
count = 0
minx = 1.0e6
maxx = -1.0e6
miny = 1.0e6
maxy = -1.0e6
minz = 1.0e6
maxz = -1.0e6

for node in mesh.iter_entities(sel, StkRank.NODE_RANK):
    count += 1
    xyz = coords.get(node)
    minx = min(minx, xyz[0])
    maxx = max(maxx, xyz[0])
    miny = min(miny, xyz[1])
    maxy = max(maxy, xyz[1])
    minz = min(minz, xyz[2])
    maxz = max(maxz, xyz[2])

print("Num. nodes = ", count, "; minx = ", minx, "; maxx = ", maxx)
print("Num. nodes = ", count, "; miny = ", miny, "; maxy = ", maxy)
print("Num. nodes = ", count, "; minz = ", minz, "; maxz = ", maxz)

## Define the 2D ridge parameters
delta = 0.26
h=0.04/delta;
b=0.1/delta;
Lx=maxx-minx;
Ly=maxy-miny;
Lz=maxz-minz;
a=3.; # Mesh damping factor
## Interating and acting over entire buckets
for bkt in mesh.iter_buckets(sel, StkRank.NODE_RANK):
    xyz = coords.bkt_view(bkt)
    mask=np.heaviside(b-np.abs(xyz[:,0]),1.0)
    xyz[:, -1] = xyz[:,-1]+mask*h*np.cos(np.pi*xyz[:,0]/(2*b))*(1-xyz[:,-1]/Lz)**a


# Iterating over entities
count = 0
miny = 1.0e6
maxy = -1.0e6
for node in mesh.iter_entities(sel, StkRank.NODE_RANK):
    count += 1
    xyz = coords.get(node)
    miny = min(miny, xyz[2])
    maxy = max(maxy, xyz[2])

print("Num. nodes = ", count, "; minz = ", minz, "; maxz = ", maxz)

filename_out="2d_ridge.exo"
purpose=stk.DatabasePurpose.WRITE_RESULTS
stkio = mesh.stkio
fh = stkio.create_output_mesh(filename_out, purpose=purpose)

stkio.begin_output_step(fh, 1)
stkio.write_defined_output_fields(fh)
stkio.end_output_step(fh)

del mesh
par.finalize()
