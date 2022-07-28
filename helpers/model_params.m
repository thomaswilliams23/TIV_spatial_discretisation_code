%% grid setup
cells_wide=120;
cells_long=120;
total_cells = cells_wide * cells_long;

periodic_BCs = 1;

%% simulation setup
dt=0.05;
final_time=2000;

MAX_TIME_UNINFECTED=10; %amount of time the system can spend without
                        %any active infected cells before we call it dead

%% initial values
moi = 0.01;

T_tot=4e8;     %this is the total number of cells in the respiratory tract, used to derive beta, etc
prop_of_tract = T_tot/total_cells;  %used to weight contact parameters for the size of this patch

%% visualisation
vis_grid=0;
vis_net=0;
vis_virus=0;

vis_dt=1/dt;

%% model parameters
beta_param = prop_of_tract*1.58e-8/24.0;%from Hernandez-Vargas and Velasco-Hernandez
p = 5.36/24.0;           %ibid
c = 2.4/24.0;            %ibid
delta_param = 1.04/24.0;   %ibid

virus_diff=0.1;            %CAUTION: if increasing, need to shrink dt

%% saving output
save_time_series=1;
save_params=1;
save_plots=1;




