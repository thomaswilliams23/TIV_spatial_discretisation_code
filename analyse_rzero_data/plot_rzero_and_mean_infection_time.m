%%% analyses data on infection times from an r_zero experiment and 
%   generates a plot of R_0, and also mean infection time for each dx 

source_folder = 'Results';

%ENSURE MATCHES WITH PARAMETERS USED
grid_size = 120^2;
moi_val = 0.01;

num_init_infected = round(moi_val*grid_size);

dx_vals = [1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6];
MAX_INFECTIONS = 20000;

cell_lifetime = 1/(1.04/24.0);
%%%




%initialise
mean_infection_time = zeros(length(dx_vals), 1);
mean_infection_sd = zeros(length(dx_vals), 1);

r_zero = zeros(length(dx_vals), 1);
r_zero_sd = zeros(length(dx_vals), 1);


%loop over dx values
for run_num = 1:length(dx_vals)
    
    %open files
    fhandle = load(strcat(source_folder, '/Sweep_run_', num2str(run_num), '/infection_times.mat'));
    infection_times = fhandle.infection_times;
    
    
    %work out overall mean infection time [MEAN OF AVERAGE INFECTION TIME FOR EACH SIM]
    mean_inf_time_estimates = zeros(1,size(infection_times,1));
    for inst = 1:size(infection_times,1)
        inf_times_this_inst = infection_times(inst,:);
        mean_inf_time_estimates(inst) = mean(inf_times_this_inst(inf_times_this_inst>0));
    end
    mean_infection_time(run_num) = mean(mean_inf_time_estimates);
    mean_infection_sd(run_num) = std(mean_inf_time_estimates);
        
    
    
    %work out rzero 
    r_zero(run_num) = sum(sum((infection_times>0)))/(size(infection_times, 1)*num_init_infected);
    r_zero_sd(run_num) = std(sum((infection_times>0),2)/num_init_infected);
 
end

%trim mean time of infection
actual_max_infections = max(sum((mean_time_of_nth_infection>0), 2));
mean_time_of_nth_infection = mean_time_of_nth_infection(:, 1:actual_max_infections);

mean_time_of_nth_infection(mean_time_of_nth_infection==0) = NaN;


%% plot mean infection time
figure
errorbar(log10(dx_vals), mean_infection_time, mean_infection_sd);
xticks(log10(dx_vals))
xticklabels({'^1/_6', '^1/_5', '^1/_4', '^1/_3', '^1/_2', '1', '2', '3', '4', '5', '6'})
xlabel('log_{10}(\Delta x)')
ylabel('mean infection time (h)')


%% plot rzero
figure
errorbar(log10(dx_vals), r_zero, r_zero_sd);
xticks(log10(dx_vals))
xticklabels({'^1/_6', '^1/_5', '^1/_4', '^1/_3', '^1/_2', '1', '2', '3', '4', '5', '6'})
xlabel('\Delta x (log scale)')
ylabel('R_0^*')
