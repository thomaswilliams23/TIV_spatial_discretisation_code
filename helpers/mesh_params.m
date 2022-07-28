
coarse_mesh = (pde_dx>1);

mesh_width = round(cells_wide/pde_dx)+2;
mesh_length = round(cells_long/pde_dx)+2;
total_mesh_nodes = mesh_width*mesh_length;

courant = virus_diff*pde_dt/pde_dx^2;
assert(courant<=0.25)