function local_virus = get_local_density(field, x, y, step_width)
%finds the amount of "field" in the vicinity of cell (x,y)



if step_width<=1 %fine or coincident mesh
    
    %local virus is average node height in the domain of the cell * cell
    %area (which is just 1)
    
    local_virus=0;
    for i=1:round(1/step_width)
        for j=1:round(1/step_width)
            local_virus = local_virus+field(round((x-1)/step_width) + 1 + i,...
                round((y-1)/step_width) + 1 + j) * step_width^2;
        end
    end
    
else %coarse mesh
    
    cell_area=1.0^2;
    local_virus = field(ceil(x/step_width)+1,ceil(y/step_width)+1) * cell_area;
end
