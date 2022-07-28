%%% uses output from an rzero simulation to make an infection over time
%%% plot as in Figure 6a



target_folder = '../Results/Sweep_run_1';  %should specify a simulation and dx value

fhandle = load(strcat(target_folder, '/infection_times.mat'));
infection_times = fhandle.infection_times;

rep_num = 1; %only do the one rep too

dt=0.05;
final_time = 80;

infections_over_time = zeros(ceil(final_time/dt) + 1, 1);

for j = 1:size(infection_times,2) %loop over potential infections
    
    %check if there was a jth infection
    if infection_times(rep_num,j)>0
        
        %count cumulative number of infections
        for time_ind = ceil(infection_times(rep_num,j)/dt):ceil(final_time/dt)+1
            infections_over_time(time_ind) = infections_over_time(time_ind) + 1;
        end
        
    end
end

%plot infections
plot(0:dt:final_time, infections_over_time)
hold on
plot([0, final_time], [infections_over_time(end), infections_over_time(end)], 'k--')

%add time of cell death
delta_param = 1.04/24.0;
mean_cell_lifetime = 1/delta_param;
plot([mean_cell_lifetime, mean_cell_lifetime], [0, 1.2*infections_over_time(end)], 'k:')

ylim([0, 1.2*infections_over_time(end)])

xticks([])
yticks([])

xlabel('time')
ylabel('num. infections')

box off



