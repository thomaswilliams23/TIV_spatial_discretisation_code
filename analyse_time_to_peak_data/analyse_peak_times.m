%%% Analyse peak viral load times from a simulation sweep and make a plot
%%% as in Figure 7c of the paper


%source folder
source_folder = '../Results';

%size params
dx_sweep = [1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6];
num_dx = length(dx_sweep);

num_reps = 200;

virus_peak_times = zeros(num_dx, num_reps);


%read data, get peak times
for dx=1:num_dx

    data_handle = open(strcat(source_folder, '/Sweep_run_', num2str(dx), '/net_V.mat'));
    virus_time_series = data_handle.net_V;
    [~, peak_times] = max(virus_time_series,[], 1);

    virus_peak_times(dx,:) = peak_times;
end



%find mean non-dieout time to peak
mean_peaks = 0*dx_sweep;
peak_std = 0*dx_sweep;
dieout_cutoff = 50;
for dx=1:num_dx
    non_dieout_ind = find(virus_peak_times(dx,:)>dieout_cutoff);
    non_dieout =  virus_peak_times(dx, non_dieout_ind);
    
    mean_peaks(dx) = 0.1*mean(non_dieout);
    peak_std(dx) = 0.1*std(non_dieout);
end



%plot
figure

errorbar(log10(dx_sweep), mean_peaks, peak_std);


xticks(log10(dx_sweep))
xticklabels({'^1/_6', '^1/_5', '^1/_4', '^1/_3', '^1/_2', '1', '2', '3', '4', '5', '6'})


xlabel('\Delta x (log scale)')
ylabel('time to peak (h)')



