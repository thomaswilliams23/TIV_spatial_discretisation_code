function prod_matrix = get_production_matrix(prod_cell, cell_grid, stepwidth, weight)
%makes a matrix the same size as (the interior of) the PDE grids, which is
%n at nodes which contain n producing cells and zero if there are none. For
%fine meshes, reports a 1 if the node corresponds to part of a producing
%cell

%initialise
prod_matrix = zeros(round(size(cell_grid,1)/stepwidth));
is_producing = (cell_grid==prod_cell);

%loop over cells and increment
for i=1:size(is_producing,1)
    for j=1:size(is_producing,2)
        if is_producing(i,j)
            
            if stepwidth<=1 %fine or coincident mesh
                for x = 1:round(1/stepwidth)
                    for y = 1:round(1/stepwidth)
                        prod_matrix(round((i-1)/stepwidth)+x, round((j-1)/stepwidth)+y) = ...
                            prod_matrix(round((i-1)/stepwidth)+x, round((j-1)/stepwidth)+y) + weight(i,j);
                    end
                end
                
            else %coarse mesh
                prod_matrix(ceil(i/stepwidth), ceil(j/stepwidth)) = ...
                    prod_matrix(ceil(i/stepwidth), ceil(j/stepwidth)) + weight(i,j);
            end
        end
    end
end
