%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       VIRUS PROBE SWEEP
%
%
% Produces Figures 4 and 5 from the paper. For small dx values, runs the
% simplified model once; for large dx values, runs the model multiple times
% and averages out results. Also keeps track of total virus in the system.
%
% Figure 3 from the paper can be reproduced as follows: 
% - modify dx_sweep to a single value (we used dx=1/6)
% - set diff_D = 10 and dt = 0.0005 (to ensure convergence)
% - fix source_c = [30,30] at line 115
% - change the "parfor" loop to a "for" loop
% - call contourf on the pde_mesh after simulation
%
% SI Figures 1 and 2 can be reproduced by modifying the diffusion
% coefficient diff_D and adjusting dt accordingly to ensure convergence


addpath helpers


%% sweep setup
dx_sweep = [1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6];
num_sweep=length(dx_sweep);


%% params
%computational
dt=0.05;%6e-4;
final_time=50;

nreps = 200; %for coarse mesh case

num_steps = round(final_time/dt)+1;

save_frame = 1;

cells_wide=120;
cells_long=120;
total_cells=cells_wide*cells_long;

centre_cell_x = ceil(cells_wide/2);
centre_cell_y = ceil(cells_long/2);

%model
p = 5.36/24.0;
c = 2.4/24.0; 
diff_D = 0.1;

delta_param = 1.04/24.0;
mean_cell_life = 1/delta_param;

%model options
cont_release = 1;  %if ON, sim starts with no virus but a source term
cut_off_after_cell_death = 1; % if ON, source cell stops producing virus after death 

%initial conditions
point_source_height = 1e3;

%graph
vis_grid=1;


%set up probes
probe_diag_distances = [0, 1, 30];
num_probes = length(probe_diag_distances);


%initialise
cell_probe_heat = zeros(length(dx_sweep), num_steps, length(probe_diag_distances));
total_virus_in_system = zeros(length(dx_sweep), num_steps);


%% sweep loop    
for sweep_ind=1:num_sweep

    %print header
    fprintf('Running dx value %d of %d\n', sweep_ind, num_sweep)

    %initialise
    dx=dx_sweep(sweep_ind);
    coarse_mesh = (dx>1);

    periodic_BCs=1;

    mesh_width = round(cells_wide/dx)+2;
    mesh_length = round(cells_long/dx)+2;
    total_mesh_nodes = mesh_width*mesh_length;

    nodes_per_cell = (round(1/dx))^2;

    %courant stability
    courant = diff_D*dt/dx^2;
    assert(courant<=0.25)
    
    
    
    %decide if fine or coarse mesh: this determines whether to do multiple
    %iterations
    if ~coarse_mesh
        nreps = 1;
    else
        nreps = 50;
    end
    
    
    %set up local result arrays
    virus_at_probes_this_dx = zeros(num_steps, num_probes, nreps);
    total_virus_this_dx = zeros(num_steps, nreps);
    
    
    
    %run simulations
    for rep=1:nreps
        
        
        %DEFINE SOURCE AND PROBES
        source_c = [randi(cells_wide), randi(cells_long)];
        probe_c = zeros(num_probes, 2);
        for probe_num = 1:num_probes
            probe_c(probe_num, :) = [mod(source_c(1)+probe_diag_distances(probe_num)-1, cells_wide)+1,...
                                     mod(source_c(2)+probe_diag_distances(probe_num)-1, cells_long)+1];
        end
        
        

        %make cell grid
        cells = zeros([cells_wide, cells_long]);
        cells(source_c(1), source_c(2))=1;


        %make PDE grid
        pde_mesh = zeros([mesh_width, mesh_length]);


        %if cont_release OFF, start the simulation with a fixed amount of virus
        if ~cont_release
            if ~coarse_mesh %fine mesh
                for i = 1:round(1/dx)
                    for j = 1:round(1/dx)
                        pde_mesh(round((source_c(1)-1)/dx)+1+i,round((source_c(2)-1)/dx)+1+j)=point_source_height;
                    end
                end
            else %coarse mesh
                cells_per_node = dx^2;
                pde_mesh(ceil(source_c(1)/dx)+1, ceil(source_c(2)/dx)+1)=point_source_height/cells_per_node;
            end
        end



        %% simulate
        for time_ind = 1:num_steps

            next_mesh = pde_mesh;

            %production
            if cont_release

                prod_matrix = get_production_matrix_at_coords(cells, source_c, dx);
                vol_weight = 1 + coarse_mesh * (-1 + dx^2);

                %case where cell has died
                if cut_off_after_cell_death && (time_ind*dt > mean_cell_life)
                    prod_matrix = 0*prod_matrix;
                end

                next_mesh(2:mesh_width-1, 2:mesh_length-1) = ...
                    next_mesh(2:mesh_width-1, 2:mesh_length-1) + ...
                    p*(dt/vol_weight)*prod_matrix;
            end


            %diffusion on interior of grid
            next_mesh(2:mesh_width-1, 2:mesh_length-1) = ...
                next_mesh(2:mesh_width-1, 2:mesh_length-1) + ...
                (diff_D*dt/dx^2) * (pde_mesh(1:mesh_width-2, 2:mesh_length-1) +...
                pde_mesh(3:mesh_width, 2:mesh_length-1) + ...
                pde_mesh(2:mesh_width-1, 1:mesh_length-2) + ...
                pde_mesh(2:mesh_width-1, 3:mesh_length) - ...
                4*pde_mesh(2:mesh_width-1, 2:mesh_length-1));


            %update the ghost nodes
            if periodic_BCs
                next_mesh(1, 2:end-1) = next_mesh(end-1,2:end-1);
                next_mesh(end, 2:end-1) = next_mesh(2, 2:end-1);
                next_mesh(2:end-1,1) = next_mesh(2:end-1,end-1);
                next_mesh(2:end-1,end) = next_mesh(2:end-1,2);
            else
                next_mesh(1, 2:end-1) = next_mesh(2,2:end-1);
                next_mesh(end, 2:end-1) = next_mesh(end-1, 2:end-1);
                next_mesh(2:end-1,1) = next_mesh(2:end-1,2);
                next_mesh(2:end-1,end) = next_mesh(2:end-1,end-1);
            end

            pde_mesh=next_mesh;


            
            
            %% record virus at probe cell
            if ~mod(round(time_ind*dt), save_frame)
                virus_at_cell = zeros(1,num_probes);

                if ~coarse_mesh   %FINE MESH VERSION
                    mesh_element_area = dx^2;

                    %loop over the cells to record at
                    for cell_ind = 1:size(probe_c, 1)
                        for i = 1:round(1/dx)
                            for j = 1:round(1/dx)
                                virus_at_cell(cell_ind) = virus_at_cell(cell_ind) + ...
                                    pde_mesh(round((probe_c(cell_ind,1)-1)/dx)+1+i,...
                                             round((probe_c(cell_ind,2)-1)/dx)+1+j)*...
                                    mesh_element_area;
                            end
                        end
                    end
                    
                else    %COARSE MESH VERSION
                    for cell_ind = 1:size(probe_c, 1)
                        virus_at_cell(cell_ind) = pde_mesh(ceil(probe_c(cell_ind,1)/dx)+1,...
                            ceil(probe_c(cell_ind,2)/dx)+1);
                    end
                end

                virus_at_probes_this_dx(time_ind, :, rep) = virus_at_cell;
                total_virus_this_dx(time_ind, rep) = sum(sum(pde_mesh(2:end-1,2:end-1)))*dx^2;

            end
        end        
    end
 
    %pass out info
    cell_probe_heat(sweep_ind, :, :) = mean(virus_at_probes_this_dx, 3);
    total_virus_in_system(sweep_ind, :) = mean(total_virus_this_dx, 2);
end



%% plot results

for cell_ind = 1:num_probes

    %time series plot:
    
    %setup
    plot_handles_time_series = [];
    plot_labels_small = {'1/6', '1/5', '1/4', '1/3', '1/2', '1'};
    plot_labels_large = {'1', '2', '3', '4', '5', '6'};

    
    %define the colours
    colour_set = jet(11);
    
    %small dx
    colour_set(1,:) = [0, 0, 0.67];
    colour_set(2,:) = [0, 0, 1];
    colour_set(3,:) = [0, 0.33, 1];
    colour_set(4,:) = [0, 0.67, 1];
    colour_set(5,:) = [0, 1, 1];
    
    %reference
    colour_set(6, :) = [0,0,0];
    
    %big dx
    colour_set(7,:) = [1, 0.85, 0];
    colour_set(8,:) = [1, 0.7, 0];
    colour_set(9,:) = [1, 0.5, 0];
    colour_set(10,:) = [1, 0.25, 0];
    colour_set(11,:) = [1, 0, 0];
    

    
    %small dx only
    figure
    for dx_ind = 1:6
        plot_handles_time_series(dx_ind) = plot(dt*(1:num_steps), squeeze(cell_probe_heat(dx_ind, :, cell_ind)),...
            'Color', colour_set(dx_ind, :));
        hold on
    end
    
    %dead cell time
    %plot([1/delta_param, 1/delta_param], [0, 0.9], 'k--');
    
    hold off

    legend(plot_handles_time_series, plot_labels_small, 'Location', 'NorthEast')
    xlim([0,dt*num_steps])
    xlabel('time (hours)')
    ylabel('virus at probe cell')
    
    
    
    
    %large dx only
    figure
    for dx_ind = 6:11
        plot_handles_time_series(dx_ind-5) = plot(dt*(1:num_steps), squeeze(cell_probe_heat(dx_ind, :, cell_ind)),...
            'Color', colour_set(dx_ind, :));
        hold on
    end
    uistack(plot_handles_time_series(1), 'top');
    
    %dead cell time
    %plot([1/delta_param, 1/delta_param], [0, 0.9], 'k--');
        
    
    hold off

    legend(plot_handles_time_series, plot_labels_large, 'Location', 'NorthEast')
    xlim([0,dt*num_steps])
    xlabel('time (hours)')
    ylabel('av. virus at probe cell')
    

end




%and prop virus exported


%small dx
figure
prop_export_handles = zeros(6, 1);
for dx_ind = 1:6
    prop_export_handles(dx_ind) = plot(dt*(1:num_steps), 1 - squeeze(cell_probe_heat(dx_ind, :, 1))./total_virus_in_system(dx_ind, :),...
        'Color', colour_set(dx_ind, :));
    hold on
end

%dead cell time
plot([1/delta_param, 1/delta_param], [0, 1], 'k--');

hold off

legend(prop_export_handles, plot_labels_small, 'Location', 'SouthEast')
xlim([0,dt*num_steps])
ylim([0,1])
xlabel('time (hours)')
ylabel('prop. virus exported')




%large dx
figure
prop_export_handles = zeros(6, 1);
for dx_ind = 6:11
    prop_export_handles(dx_ind-5) = plot(dt*(1:num_steps), 1 - squeeze(cell_probe_heat(dx_ind, :, 1))./total_virus_in_system(dx_ind, :),...
        'Color', colour_set(dx_ind, :));
    hold on
end
uistack(prop_export_handles(1), 'top');

%dead cell time
plot([1/delta_param, 1/delta_param], [0, 1], 'k--');

hold off

legend(prop_export_handles, plot_labels_large, 'Location', 'SouthEast')
xlim([0,dt*num_steps])
ylim([0,1])
xlabel('time (hours)')
ylabel('av. prop. virus exported')




    
%also plot total virus
figure
plot(dt*(1:num_steps), total_virus_in_system(dx_ind, :), 'k');
xlim([0,dt*num_steps])
xlabel('time (hours)')
ylabel('total virus')

