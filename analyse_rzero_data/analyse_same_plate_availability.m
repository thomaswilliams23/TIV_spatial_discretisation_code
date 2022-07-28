%%% Analyses the spatial data from an rzero simulation and generates a
%%% figure as in SI Figure 4d which charts the proportion of cells sharing
%%% a node with an infected cell which become infected over time



folder_stem = '../Results';

target_folders = [2,6]; %specify which dx values to analyse here
dx_vals = target_folders;


%experiment specs   %%%%%%%%%%%%%%%%
final_time = 80;
dt = 0.5;
moi=0.01;
grid_size = 120^2;
num_init_infected = moi*grid_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for folder_num = 1:length(target_folders)
    
    %initialise time series
    num_same_plate_infected = zeros(round(final_time/dt),1);
    
    
    
    
    figure
    ax=axes;
    
    %open data
    fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folders(folder_num)), '/infection_times.mat'));
    infection_times = fhandle.infection_times;
    
    fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folders(folder_num)), '/dist_from_infected_plate.mat'));
    plate_distances = fhandle.dist_from_infected_plate;
    
    fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folders(folder_num)), '/num_same_plate_cells.mat'));
    num_same_plate_cells = fhandle.num_same_plate_cells;
    
    
    %initialise vectors for infection times
    inf_times_as_vec = zeros(sum(sum((infection_times>0))),1);
    same_plate_inf_times_as_vec = 0*inf_times_as_vec;
    
    
    %loop over all infection events
    v_1 = 1;
    v_2 = 1;
    for i = 1:size(infection_times,1) %reps

        for j = 1:size(infection_times,2) %infections
            if infection_times(i,j)>0
                
                %save all infection times into a vector
                inf_times_as_vec(v_1) = infection_times(i,j);
                
                %also set aside same-plate infections
                if plate_distances(i,j)==0
                    same_plate_inf_times_as_vec(v_2) = infection_times(i,j);
                    
                    %count cumulative number of infections
                    for time_ind = ceil(infection_times(i,j)/dt):ceil(final_time/dt)
                        num_same_plate_infected(time_ind) = ...
                            num_same_plate_infected(time_ind) + 1;
                    end
                    
                    v_2 = v_2 + 1;
                end
                
                v_1 = v_1 + 1;
            end
        end
    end
    same_plate_inf_times_as_vec = same_plate_inf_times_as_vec(2:v_2-1); %trim extra length
    
    
    %rescale time series
    prop_same_plate_infected = num_same_plate_infected/sum(num_same_plate_cells);
    
    
    %plot
    plot(dt:dt:final_time, prop_same_plate_infected)
    xlabel('time (h)')
    ylabel('prop. same-plate cells infected')
    
    ylim([0,1])
    xlim([0,50])
    
    hold on
    
    %also plot time of cell death     
    plot(23.0769+0*linspace(0,1), linspace(0,1), 'k:');
    
    
    
    %now calculate the point of inflection (gives an idea of rate of
    %infection decreasing)
    prop_diff = diff(prop_same_plate_infected);
    [~,max_ind] = max(smoothdata(prop_diff, 'gaussian', 10));
    
    
    POI_colour = [1, 0, 1];%(1/256)*[256, 173, 32];
    plot(max_ind*dt, prop_same_plate_infected(max_ind), 'x', 'Color', POI_colour)
    
    plot([max_ind*dt, max_ind*dt], [0, prop_same_plate_infected(max_ind)], '--', 'Color', POI_colour)
    
    %text
    tb = annotation('textbox', [0.5, 0.8, 0.1, 0.1], 'String', '1/\delta');
    tb.FontSize = 9;
    tb.EdgeColor = 'none';
    
end

