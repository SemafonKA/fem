#pragma once

/// NOTE: This file includes ALL FEM-related includes at once
/// Please INCLUDE THIS file and then use namespaces (like `using namespace fem::two_dim`) 
/// to select what realization you need.
/// In general, you can include only few of headers below, but for what?


/// DOMAINS

#include "../two_dim/two_dim_domain.h"
#include "../three_dim/three_dim_domain.h"



/// POINTS
#include "../two_dim/two_dim_point.h"
#include "../three_dim/three_dim_point.h"



/// MESHES
#include "../two_dim/two_dim_quads_linear_mesh.h"
#include "../three_dim/three_dim_cuboid_linear_mesh.h"



/// GRIDS
#include "../two_dim/two_dim_quads_linear_grid.h"
#include "../three_dim/three_dim_cuboid_linear_grid.h"



/// SOLVERS
#include "../two_dim/two_dim_quads_linear_solver.h"
#include "../three_dim/three_dim_cuboid_linear_solver.h"
