%%% Analyse the output from a model simulation sweep and plot multiple time
%%% series for net target and infected cells and virus as in Figure 1 of
%%% the paper




%experiment specs
dt=0.05;
save_dt = 0.1/dt;
%


%number reps to plot
num_time_series = 10;


%open data
source_folder = '../Results/Sweep_run_1'; %need to specify a folder and dx value

tmp_handle = open(strcat(source_folder, '/net_T'));
target_time_series = tmp_handle.net_T;

tmp_handle = open(strcat(source_folder, '/net_I'));
infected_time_series = tmp_handle.net_I;

tmp_handle = open(strcat(source_folder, '/net_V'));
virus_time_series = tmp_handle.net_V;




%plot everything
figure

hold on


time_vec=((dt*save_dt)*(1:size(net_V,1)))';

%colours
target_colour = [0 0.447 0.741];
infected_colour = [0.85 0.325 0.098];
virus_colour = [0.9290 0.6940 0.1250];

%plot each time series
for i = 1:num_time_series
    plot(time_vec,net_T(:,i),'Color',target_colour, 'LineWidth', 0.1, 'LineStyle', '-');
    plot(time_vec,net_I(:,i),'Color',infected_colour, 'LineWidth', 0.1, 'LineStyle', '-');
    plot(time_vec,net_V(:,i),'Color',virus_colour, 'LineWidth', 0.1, 'LineStyle', '-');
end


%plot mean
plot_t=plot(time_vec,mean(net_T(:,1:num_time_series),2),'Color',target_colour, 'LineWidth', 2);
plot_i=plot(time_vec,mean(net_I(:,1:num_time_series),2),'Color',infected_colour, 'LineWidth', 2);
plot_v=plot(time_vec,mean(net_V(:,1:num_time_series),2),'Color',virus_colour, 'LineWidth', 2);


xlim([0,500])
ylim([0,1.1])

legend([plot_t, plot_i, plot_v], {'T/N', 'I/N', 'V/V_{max}'})

xlabel('hours post infection')
ylabel('cell fraction / V/V_{max}')

hold off
