%%% Analyses the spatial data from an rzero simulation and generates a
%%% figure as in SI Figure 4a which plots a histogram of the time of
%%% infections in a given simulation, stratified by the distance between
%%% the newly infected cell and the nearest already-infected cell.



folder_stem = '../Results';
target_folder = 1; %choice of dx value, dx=1 works best


figure
ax=axes;

%open data
fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folder), '/infection_times.mat'));
infection_times = fhandle.infection_times;

fhandle = load(strcat(folder_stem, '/Sweep_run_', num2str(target_folder), '/dist_from_infected_plate.mat'));
plate_distances = fhandle.dist_from_infected_plate;


%initialise vectors for infection times
inf_times_as_vec = zeros(sum(sum((infection_times>0))),1);

%make some arrays for infection times within certain radii
radius_vals = 1:6;
for radius_ind = 1:length(radius_vals)
    at_least_r_from_source{radius_ind} = 0*inf_times_as_vec;
end


%loop over all infection events
vec_ind = 1;
thresh_vec_inds = 1+0*radius_vals;
for i = 1:size(infection_times,1)
    for j = 1:size(infection_times,2)
        if infection_times(i,j)>0

            %save all infection times into a vector
            inf_times_as_vec(vec_ind) = infection_times(i,j);

            %also set aside infections within each radius
            for radius_ind = 1:length(radius_vals)
                if plate_distances(i,j)<=radius_vals(radius_ind)
                    at_least_r_from_source{radius_ind}(thresh_vec_inds(radius_ind)) = infection_times(i,j);
                    thresh_vec_inds(radius_ind) = thresh_vec_inds(radius_ind) + 1;
                end
            end

            vec_ind = vec_ind + 1;
        end
    end
end

%plot as histogram
hist_handles = zeros(length(radius_vals)+1,1);
hist_labels{1} = 'all';

hist_handles(1) = histogram(ax, inf_times_as_vec, 'BinWidth', 1);

max_height = 70*size(infection_times,1);
ylim([0,max_height])

xlabel('time of infection (h)')
ylabel('freq.')

hold on

for radius_ind = length(radius_vals):-1:1
    same_plate_inf_times_as_vec = at_least_r_from_source{radius_ind}(1:thresh_vec_inds(radius_ind)-1);
    
    hist_handles(2+length(radius_vals)-radius_ind) = histogram(ax, same_plate_inf_times_as_vec, 'BinWidth', 1);
    
    %make a label too
    hist_labels{2+length(radius_vals)-radius_ind} = sprintf('within %d CD', radius_vals(radius_ind));
end

%also plot time of cell death
plot([23.0769, 23.0769], [0,max_height], 'k:');

%also plot mean time of infection
mean_inf_time = mean(inf_times_as_vec);
plot([mean_inf_time, mean_inf_time], [0,max_height], 'r--')

legend(hist_handles, hist_labels)



