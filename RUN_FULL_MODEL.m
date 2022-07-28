%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  TIV CELLULAR AUTOMATA SIMULATION
%
% 
% This code implements the full multiscale version of the TIV model as
% described in the paper. The model is called nreps times in parallel for
% each of the dx values specified in dx_sweep.
%
% This code in its current form is set up to perform time-to-peak 
% simulations as in the Results section of the paper. When called, it saves
% its output to a "Results" folder. This output can be analysed by the code 
% in folder "analyse_time_to_peak_data", which generates the plots in 
% Figure 7c of the paper. To generate graphics of the state of the cell
% grid as in Figure 1 and 7a of the paper, change the "parfor" loop to a 
% for loop, then specify the output times under the OPTIONS heading below
% and uncomment the code section at line 236.
%
% This code is easily modified to reproduce Figure 1 in the paper by making
% the following parameter reassingnments:
% - moi = 1/total_cells;
% - final_time = 700; (for virus_diff = 0.1, readjust for other diffusion)
% - nreps = 10;
% - change virus_diff and dt as required
% Then specify grid_vis_times as instructed above and call the code in the
% folder "analyse_model_demo_output" on the model output.
% The folder "analyse_model_demo_output" also contains the code for the ODE
% form of the model and generates the bottom row of graphics for Figure 1.




addpath helpers

%initialise
model_params;

nreps = 200;
final_time = 250; %sufficient to obtain time of peak

save_dt = 0.1/dt;


%%%% OPTIONS

%to visualise the state of the cell grid at specified times, uncomment
%below
%grid_vis_times = [0, 30, 60, 90];

%default diffusion is 0.1, to change to 1 uncomment below
%virus_diff = 1;
%dt = 0.005; %must satisfy Courant stability

%default moi is 0.01, to change to a single initially-infected cell (as in 
%Figure 1) uncomment below
%moi=1/total_cells;

%%%%



%initialise net dynamics arrays
net_T = zeros([floor((round(final_time/dt)+1)/save_dt)+1,nreps]);
net_I = zeros([floor((round(final_time/dt)+1)/save_dt)+1,nreps]);
net_V = zeros([floor((round(final_time/dt)+1)/save_dt)+1,nreps]);


%dx values
dx_sweep = [1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6];
num_dx=length(dx_sweep);



%set up results folder
mkdir('Results');


fprintf('\n\n____________RUNNING PARAMETER SWEEP____________\n');
%sweep over parameters
for sweep_ind=1:num_dx

    fprintf('\n SIMULATIONS WITH dx = %.3f\n',dx_sweep(sweep_ind));

    %params
    dx=dx_sweep(sweep_ind);
    mesh_params;

    %set up output folder
    result_folder = sprintf('Sweep_run_%d',sweep_ind);
    mkdir('./Results/',result_folder);
    dest_folder = strcat('./Results/', result_folder, '/');
    
    

    %% run sim nreps times in parallel
    parfor rep=1:nreps

        %% initialise
        %virus
        mesh_weight = 1 + (dx>1)*(dx^2 - 1); %weight to account for virus having 
                                             %to spread across a wide element
        virus = zeros([mesh_width,mesh_length]);

        %initialise the cell grid
        grid=zeros([cells_wide, cells_long]);

        %now seed the grid with some number of infected cells in random
        %locations
        for inf_cell = 1:round(moi*total_cells)
            x_c=randi(cells_wide);
            y_c=randi(cells_long);

            %check not already infected
            while grid(x_c,y_c)==2
                x_c=randi(cells_wide);
                y_c=randi(cells_long);
            end

            grid(x_c,y_c)=2;
        end
        
        
        %initialise time series
        T_time_series = zeros([1,round(final_time/dt)+1]);
        I_time_series = zeros([1,round(final_time/dt)+1]);
        D_time_series = zeros([1,round(final_time/dt)+1]);
        V_time_series = zeros([1,round(final_time/dt)+1]);

        %put in starting value
        T_time_series(1)=sum(sum((grid==1)))/total_cells;
        I_time_series(1)=sum(sum((grid==2)))/total_cells;
        D_time_series(1)=sum(sum((grid==0)))/total_cells;
        V_time_series(1)=sum(sum((virus)));





        %% main loop
        t=0;
        time_last_infected=0;
        infection_complete=0;
        infection_dead=0;
        while t<=final_time && ~infection_dead && ~infection_complete

            %% update grid
            %loop over cells
            new_grid=grid;
            for i=1:cells_wide
                for j=1:cells_long

                    %check cell type:
                    if grid(i,j)==1  %cell is TARGET

                        %compute probability of being infected on next time step
                        prob_t_to_i = 1-exp(-total_cells*beta_param*get_local_density(virus,i,j,dx)*dt);

                        %draw random number, decide next state of cell
                        if rand<prob_t_to_i  %becomes infected
                            new_grid(i,j)=2;
                        end
                        

                    elseif grid(i,j)==2   %cell is INFECTED

                        %compute prob of cell dying in next time step
                        prob_i_to_d = 1-exp(-delta_param*dt);

                        %draw random number, decide next state of cell
                        if rand<prob_i_to_d %becomes dead
                            new_grid(i,j)=0;
                        end


                    end  %dead cells do nothing

                end
            end              



            %% update PDE

            %update virus over cells
            new_virus=virus;

            %production volume weighting for large elements
            vol_weight = 1 + coarse_mesh * (-1 + dx^2);

            %decay
            new_virus = new_virus - dt*c*virus;

            %diff
            new_virus(2:mesh_width-1, 2:mesh_length-1) = ...
                new_virus(2:mesh_width-1, 2:mesh_length-1) + ...
                virus_diff*(dt/dx^2) * (virus(1:mesh_width-2, 2:mesh_length-1) +...
                virus(3:mesh_width, 2:mesh_length-1) + ...
                virus(2:mesh_width-1, 1:mesh_length-2) + ...
                virus(2:mesh_width-1, 3:mesh_length) - ...
                4*virus(2:mesh_width-1, 2:mesh_length-1));

            %production
            infected_dist = get_production_matrix(2, grid, dx, 0*grid+1);

            new_virus(2:mesh_width-1, 2:mesh_length-1)=...
                new_virus(2:mesh_width-1, 2:mesh_length-1)...
                + (dt/vol_weight)*p*infected_dist;

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


            %extract proportion of each cell type for time series
            T_time_series(round(t/dt)+1) = sum(sum((grid==1)))/total_cells;
            I_time_series(round(t/dt)+1) = sum(sum((grid==2)))/total_cells;
            D_time_series(round(t/dt)+1) = sum(sum((grid==0)))/total_cells;
            V_time_series(round(t/dt)+1) = sum(sum(virus(2:end-1,2:end-1)))*dx^2;
            
            
            
            %% plot state of grid [ONLY if not in parallel]
%             for vis_ind = 1:length(grid_vis_times)
%                 if abs(t-grid_vis_times(vis_ind))<0.5*dt
%                     pause(1);
%                     plot_grid;
%                     saveas(gcf, strcat('Results/grid_frame_', num2str(vis_ind)), 'fig');
%                 end
%             end



            %% check if there is any reason to end the loop
            %check if there are any infected cells remaining, if so update
            %the clock for this
            if I_time_series(round(t/dt)+1)>0
                time_last_infected=t;
            end

            %check if the tissue is completely invaded
            infection_complete = ((T_time_series(round(t/dt)+1)==0) &&...
               (I_time_series(round(t/dt)+1)==0) && (V_time_series(round(t/dt)+1)<10));

            %or has died out
            infection_dead = ((t-time_last_infected>=MAX_TIME_UNINFECTED) &&...
                (V_time_series(round(t/dt)+1)<1));

            t=t+dt;
        end

        %fill in the rest of the time series
        if t<final_time
            last_t=t-dt;
            while t<=final_time
                T_time_series(round(t/dt)+1) = T_time_series(round(last_t/dt)+1);
                I_time_series(round(t/dt)+1) = I_time_series(round(last_t/dt)+1);
                V_time_series(round(t/dt)+1) = V_time_series(round(last_t/dt)+1);
                t=t+dt;
            end
        end


        %from the time_series arrays, extract values at an interval of 
        %t=save_dt and pass these into the output arrays
        T_at_output_values=zeros(1,floor((round(final_time/dt)+1)/save_dt)+1);
        I_at_output_values=zeros(1,floor((round(final_time/dt)+1)/save_dt)+1);
        V_at_output_values=zeros(1,floor((round(final_time/dt)+1)/save_dt)+1);
        for time_index=1:floor((round(final_time/dt)+1)/save_dt)+1
            T_at_output_values(time_index) = T_time_series(round((time_index-1)*save_dt)+1);
            I_at_output_values(time_index) = I_time_series(round((time_index-1)*save_dt)+1);
            V_at_output_values(time_index) = V_time_series(round((time_index-1)*save_dt)+1)/max(V_time_series);
        end

        net_T(:,rep) = T_at_output_values;
        net_I(:,rep) = I_at_output_values;
        net_V(:,rep) = V_at_output_values;


    end
    
    %write to file
    if save_time_series
        save(strcat(dest_folder, 'net_T'),'net_T');
        save(strcat(dest_folder, 'net_I'),'net_I');
        save(strcat(dest_folder, 'net_V'),'net_V');
    end

    %format
    fprintf('\n________________________________________\n');


end
