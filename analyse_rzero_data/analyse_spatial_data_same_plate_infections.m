%%% Analyses the spatial data from an rzero simulation and generates a
%%% figure as in SI Figure 4c which plots a histogram of the time of
%%% infections in a given simulation, stratified by whether the infection
%%% was of a cell sharing a node with an infected cell or not.


folder_stem = '../Results';

target_folders = [2,6]; %choice of dx value 

for folder_num = 1:length(target_folders)
    
    figure
    ax=axes;
    
    %open data
    fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folders(folder_num)), '/infection_times.mat'));
    infection_times = fhandle.infection_times;
    
    fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folders(folder_num)), '/dist_from_infected_plate.mat'));
    plate_distances = fhandle.dist_from_infected_plate;
    
    
    %initialise vectors for infection times
    inf_times_as_vec = zeros(sum(sum((infection_times>0))),1);
    same_plate_inf_times_as_vec = 0*inf_times_as_vec;
    
    all_distances_as_vec = 0*inf_times_as_vec;
    
    
    %loop over all infection events
    v_1 = 1;
    v_2 = 1;
    for i = 1:size(infection_times,1)
        for j = 1:size(infection_times,2)
            if infection_times(i,j)>0
                
                %save all infection times into a vector
                inf_times_as_vec(v_1) = infection_times(i,j);
                
                all_distances_as_vec(v_1) = plate_distances(i,j);
                
                %also set aside same-plate infections
                if plate_distances(i,j)==0
                    same_plate_inf_times_as_vec(v_2) = infection_times(i,j);
                    v_2 = v_2 + 1;
                end
                
                v_1 = v_1 + 1;
            end
        end
    end
    same_plate_inf_times_as_vec = same_plate_inf_times_as_vec(2:v_2-1); %trim extra length
    
    
    %plot as histogram
    h_all = histogram(ax, inf_times_as_vec, 'BinWidth', 1);
    max_height = 70*size(infection_times,1);
    
    ylim([0,max_height])
    xlabel('time of infection (h)')
    ylabel('freq.')
    
    hold on
    h_same_plate = histogram(ax, same_plate_inf_times_as_vec, 'BinWidth', 1);
    
    
    %also plot time of cell death
    plot([23.0769, 23.0769], [0,max_height], 'k:');

    %also plot mean time of infection
    mean_inf_time = mean(inf_times_as_vec);
    plot([mean_inf_time, mean_inf_time], [0,max_height], 'r--')
    
    
    legend([h_all, h_same_plate], 'all', 'same plate')
end

