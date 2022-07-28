function [node_list, num_nodes] = get_nodes_for_cell(indices, stepsize)
%produces an array of node indices which correspond to cell specified by
%'indices', assumes a mesh with ghost nodes around the boundaries

if stepsize>=1 %(coarse mesh)
    num_nodes = 1;
    node_list = [ceil(indices(1)/stepsize)+1; ceil(indices(2)/stepsize)+1];
    
else %(fine mesh)
    num_nodes = round(1/stepsize)^2;
    node_list = zeros(2,num_nodes);
    
    node_count = 1;
    for x = 1:round(1/stepsize)
        for y = 1:round(1/stepsize)
            node_list(1,node_count) = round((indices(1) - 1)/stepsize) + 1 + x;
            node_list(2,node_count) = round((indices(2) - 1)/stepsize) + 1 + y;
            node_count = node_count+1;
        end
    end
    
end
            