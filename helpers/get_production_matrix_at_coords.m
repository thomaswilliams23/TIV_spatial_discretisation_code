function prod_matrix = get_production_matrix_at_coords(cell_grid, coords, stepwidth)
%makes a matrix the same size as (the interior of) the PDE grids, which is
%1 at the cells specified by "coords" and 0 otherwise.

%get relevant nodes
[node_list, num_nodes] = get_nodes_for_cell(coords, stepwidth);
node_list = node_list - 1; %subtract 1 from everything to account for lack of boundary nodes


%fill into blank matrix
prod_matrix = zeros(round(size(cell_grid,1)/stepwidth));
for node=1:num_nodes
    prod_matrix(node_list(1,node), node_list(2,node)) = 1;
end
    
end
