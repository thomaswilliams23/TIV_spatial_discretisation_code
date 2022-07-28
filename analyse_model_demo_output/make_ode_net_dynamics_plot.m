%%% Solve the ODE form of the TIV model and plot the net dynamics as in the
%%% bottom row of Figure 1 of the paper


%% parameters
%these values are taken from Hernandez-Vargas and Velasco-Hernandez (2020)
beta_prm=1.58e-8/24.0;     %rate of infection (copies/ml per hour)
delta_prm=1.04/24.0;    %death of infected cells (per hour)
c=2.4/24.0;         %viral clearance (per hour)
p=5.36/24.0; %viral production

%% initial conditions
total_cells = 4e8;

init_moi = 1/(120)^2;

T_0 = (1-init_moi)*total_cells;
I_0 = init_moi*total_cells;
V_0 = 0;

Y_0=[T_0; I_0; V_0];


%% computational parameters
dt=0.1;
final_time=200;

t_span=linspace(0,final_time,final_time+1);


%% solve the ODE
TIV_system = @(t, y) [-beta_prm*y(1)*y(3);...
                      beta_prm*y(1)*y(3)-delta_prm*y(2);...
                      p*y(2)-c*y(3)];
[t_out, y_out] = ode23(TIV_system,t_span,Y_0);



%% plot
figure
hold on

plot(t_out,y_out(:,1)/total_cells,'LineWidth', 2);
plot(t_out,y_out(:,2)/total_cells,'LineWidth', 2);
plot(t_out,y_out(:,3)/max(y_out(:,3)),'LineWidth', 2);

xlabel('hours post infection')
ylabel('cell fraction / V/V_{max}')

legend({'T/N', 'I/N', 'V/V_{max}'})


%% produce tiles coloured by proportion of cell types
target_cells = y_out(:,1)/total_cells;
infected_cells = y_out(:,2)/total_cells;

%get amounts of different cells at an early point
time1 = 31;
target1 = target_cells(time1);
infected1 = infected_cells(time1);

%and at peak infected cells
[~, time2] = max(infected_cells);
target2 = target_cells(time2);
infected2 = infected_cells(time2);



%colours
dead_colour = [0,0,0];         %dead cells are black
target_colour = [0,0,.8];      %target cells are blue
infected_colour = [1, 0, 0];   %infected cells are red


%make two singleton heatmaps with the correct colour
figure
h1=heatmap([1]);
colour1 = target1*target_colour + infected1*infected_colour + (1-target1-infected1)*dead_colour;
colormap(colour1)
h1.CellLabelColor = 'none';


figure
h2=heatmap([1]);
colour2 = target2*target_colour + infected2*infected_colour + (1-target2-infected2)*dead_colour;
colormap(colour2)
h2.CellLabelColor = 'none';





