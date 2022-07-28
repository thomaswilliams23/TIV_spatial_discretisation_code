
coarse_mesh = (dx>1);

mesh_width = round(cells_wide/dx)+2;
mesh_length = round(cells_long/dx)+2;
total_mesh_nodes = mesh_width*mesh_length;

courant = virus_diff*dt/dx^2;
assert(courant<=0.25)