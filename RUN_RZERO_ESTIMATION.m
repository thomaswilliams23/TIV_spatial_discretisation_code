%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           RZERO ESTIMATION
%
%
% This code implements a simplified version of the multicellular version of
% the TIV model where newly-infected cells are marked for infection but do
% not produce virus (as described in the paper). The model is called nreps
% times in parallel. To visualise, set nreps=1 and replace parfor with for,
% then call plot_grid in the "helpers" folder to plot the state of the cell
% grid.
% 
% This code saves its output to a "Results" folder. This output can be 
% analysed by the code in folder "analyse_rzero_data", which generates the 
% plots in Figure 6 of the paper and SI Figures 3 and 4 in the 
% supplementary material.



addpath helpers



%parameters
model_params;

nreps = 200;
save_dt = 0.1/dt;

final_time = 80;

save_spatial_info = 1;



%%%% OPTIONS

%default diffusion is 0.1, to change to 1 uncomment below
%virus_diff = 1;
%dt = 0.005;

%default grid size is 120x120, to change to 240x240 uncomment below
%cells_wide = 240;
%cells_long = 240;
%total_cells = cells_wide*cells_long;
%beta_param = beta_param*120^2/total_cells;

%%%%


%set up some sweep values
dx_sweep = [1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6];
num_dx=length(dx_sweep);


%initialise result arrays
mean_rzero = zeros(num_dx, 1);
sd_rzero = zeros(num_dx, 1);


%set up results folder
mkdir('Results');




fprintf('\n\n____________RUNNING PARAMETER SWEEP____________\n');
%sweep over parameters
for sweep_ind=1:num_dx
    
    fprintf('\n SIMULATIONS WITH dx = %.3f\n',dx_sweep(sweep_ind));

    %params
    dx=dx_sweep(sweep_ind);
    mesh_params;


    %time of cell infection
    MAX_INFECTED = 20000;
    infection_times = zeros(nreps, MAX_INFECTED);
    
    
    %location of infection relative to nearest infected cell (only for
    %coarse mesh)
    if dx>=1 && save_spatial_info
        dist_from_infected_plate = zeros(nreps, MAX_INFECTED);
        num_same_plate_cells = zeros(nreps,1);
    end


    %% run sim nreps times in parallel
    parfor rep=1:nreps

        %% initialise
        %virus
        mesh_weight = 1 + (dx>1)*(dx^2 - 1); %weight to account for virus having 
                                                     %to spread across a wide element
        virus = zeros(mesh_width, mesh_length);
        
        
        %make the grid template for virus production (only the initially infected
        %cell ever produces virus)
        grid = ones(cells_wide, cells_long);

        %now seed the grid with some number of infected cells in random
        %locations
        num_inf_cells = round(moi*total_cells);
        coords_of_inf_cell_nodes = zeros(num_inf_cells,2);
        for inf_cell = 1:num_inf_cells
            x_c=randi(cells_wide);
            y_c=randi(cells_long);

            %check not already infected
            while grid(x_c,y_c)==2
                x_c=randi(cells_wide);
                y_c=randi(cells_long);
            end

            grid(x_c,y_c)=2;
            
            
            %compute the coordinates of the corresponding viral node
            if dx>=1 && save_spatial_info
                coords_of_inf_cell_nodes(inf_cell,:) = [ceil(x_c/dx)*dx-0.5*dx, ...
                                                        ceil(y_c/dx)*dx-0.5*dx];
            end
        end
        
        
        infected_dist = get_production_matrix(2, grid, dx, 0*grid+1);
        
        
        
        
        %also calculate how many target cells are at infected nodes
        if dx>=1 && save_spatial_info
            
            [~,unique_inf_nodes,~] = unique(coords_of_inf_cell_nodes, 'rows');
            for inf_node = 1:length(unique_inf_nodes) %for each infected node
                num_inf_cells_at_node = 0;

                %find the number of infected cells at that node
                for all_nodes = 1:size(coords_of_inf_cell_nodes,1)
                    if coords_of_inf_cell_nodes(all_nodes,:)==coords_of_inf_cell_nodes(unique_inf_nodes(inf_node),:)
                        num_inf_cells_at_node = num_inf_cells_at_node + 1;
                    end
                end

                %record number of same-plate cells available
                num_same_plate_cells(rep) = num_same_plate_cells(rep) + dx^2 - num_inf_cells_at_node;
            end
        end

        
        

        %% main loop
        t=0;
        
        %initialise output arrays
        infection_number = 1;
        infection_times_this_rep = zeros(MAX_INFECTED,1);
        if dx>=1 && save_spatial_info
            dist_from_infected_plate_this_rep = zeros(MAX_INFECTED,1);
        end
        
        %time loop
        while t<=final_time

            %% update system
            %loop over cells
            new_grid=grid;
            for i=1:cells_wide
                for j=1:cells_long

                    %check we are at a target cell
                    if grid(i,j)==1

                        %compute probability of being infected on next time step
                        prob_t_to_i = 1-exp(-total_cells*beta_param*...
                            get_local_density(virus,i,j,dx)*dt);

                        %draw random number, decide next state of cell
                        if rand<prob_t_to_i  %becomes infected
                            
                            new_grid(i,j)=4; % this is just a flag for plotting purposes
                            
                            %record infection time
                            infection_times_this_rep(infection_number) = t;
                            
                            %record distance (as euclidean distance between
                            %corresponding viral nodes) to nearest infected
                            %cell
                            if dx>=1 && save_spatial_info
                                
                                %get coordinates of plate containing new
                                %infected cell
                                new_inf_node_coords = [ceil(i/dx)*dx-0.5*dx, ceil(j/dx)*dx-0.5*dx];
                                
                                %compute euclidean distance to each node
                                %containing an infected cell
                                eu_dists = zeros(size(coords_of_inf_cell_nodes,1),1);
                                for inf_cell = 1:size(coords_of_inf_cell_nodes,1)
                                    eu_dists(inf_cell) = norm(new_inf_node_coords - coords_of_inf_cell_nodes(inf_cell,:));
                                end
                                
                                %save minimum value
                                dist_from_infected_plate_this_rep(infection_number) = min(eu_dists);
                            end
                                
                                
                                
                            infection_number = infection_number + 1;
                        end
                    end

                end
            end              



            %% update PDE
            %update virus over cells
            new_virus=virus;

            %production volume weighting for large elements
            vol_weight = 1 + coarse_mesh * (-1 + dx^2);

            %decay
            new_virus = new_virus - dt*c*virus;

            %diffusion
            new_virus(2:mesh_width-1, 2:mesh_length-1) = ...
                new_virus(2:mesh_width-1, 2:mesh_length-1) + ...
                virus_diff*(dt/dx^2) * (virus(1:mesh_width-2, 2:mesh_length-1) +...
                virus(3:mesh_width, 2:mesh_length-1) + ...
                virus(2:mesh_width-1, 1:mesh_length-2) + ...
                virus(2:mesh_width-1, 3:mesh_length) - ...
                4*virus(2:mesh_width-1, 2:mesh_length-1));

            %production (ONLY HAPPENS WHEN T < AV. INF. CELL LIFESPAN)
            if t<1/delta_param
                new_virus(2:mesh_width-1, 2:mesh_length-1)=...
                    new_virus(2:mesh_width-1, 2:mesh_length-1) + (dt/vol_weight)*p*infected_dist;
            end

            %boundaries
            if ~periodic_BCs
                new_virus(1, 2:mesh_length-1) = new_virus(2,2:mesh_length-1);
                new_virus(mesh_width, 2:mesh_length-1) = new_virus(mesh_width-1, 2:mesh_length-1);
                new_virus(2:mesh_width-1,1) = new_virus(2:mesh_width-1,2);
                new_virus(2:mesh_width-1,mesh_length) = new_virus(2:mesh_width-1,mesh_length-1);
            else
                new_virus(1, 2:mesh_length-1) = new_virus(mesh_width-1,2:mesh_length-1);
                new_virus(mesh_width, 2:mesh_length-1) = new_virus(2, 2:mesh_length-1);
                new_virus(2:mesh_width-1,1) = new_virus(2:mesh_width-1,mesh_length-1);
                new_virus(2:mesh_width-1,mesh_length) = new_virus(2:mesh_width-1,2);
            end


            %update
            virus=new_virus;


            %now update the grid
            grid=new_grid;


            t=t+dt;
        end
        
        %save infection times
        infection_times(rep,:) = infection_times_this_rep;
        
        %and spatial data
        if dx>=1 && save_spatial_info
            dist_from_infected_plate_this_rep(infection_number:end) = NaN;
            dist_from_infected_plate(rep,:) = dist_from_infected_plate_this_rep;
        end

    end
    
    %mean and standard deviation of rzero
    inst_rzero = sum((infection_times>0), 2)/(moi*total_cells);
    mean_rzero(sweep_ind) = mean(inst_rzero);
    sd_rzero(sweep_ind) = std(inst_rzero);
    
    
    %print results
    fprintf('\n        ----------------\n');
    fprintf('Mean R_zero for dx=%.3f is %.3f\n', dx, mean_rzero(sweep_ind));
    fprintf('\n________________________________________\n');


    
    %save results
    dest_folder = sprintf('./Results/Sweep_run_%d', sweep_ind);
    mkdir(dest_folder);
    
    save(strcat(dest_folder, '/infection_times'), 'infection_times');
    if coarse_mesh && save_spatial_info
        save(strcat(dest_folder, '/dist_from_infected_plate'), 'dist_from_infected_plate');
        save(strcat(dest_folder, '/num_same_plate_cells'), 'num_same_plate_cells');
    end

end


%to generate plots, call analyse_and_get_rzero_plots on the results folder




