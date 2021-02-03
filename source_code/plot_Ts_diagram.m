function [] = plot_Ts_diagram(cycle_data,my_filename,my_path,save)
%% Load the cycle points
s_cycle = [cycle_data.cycle_states.s]/1000;
T_cycle = [cycle_data.cycle_states.T]-273.15;
working_fluid = cycle_data.fluids.working_fluid;
T_min = cycle_data.properties.T_min-273.15;
T_trip = cycle_data.properties.T_trip-273.15;
T_crit = cycle_data.properties.T_crit-273.15;


%% Load the trajectories in the components
% Evaporator
s_evap_c = cycle_data.evaporator.s_c/1000;
T_evap_c = cycle_data.evaporator.T_c-273.15;
T_evap_h = cycle_data.evaporator.T_h-273.15;

% Condenser
s_cond_h = cycle_data.condenser.s_h/1000;
T_cond_h = cycle_data.condenser.T_h-273.15;
T_cond_c = cycle_data.condenser.T_c-273.15;

% Recuperator
s_rec_c = cycle_data.recuperator.s_c/1000;
s_rec_h = cycle_data.recuperator.s_h/1000;
T_rec_c = cycle_data.recuperator.T_c-273.15;
T_rec_h = cycle_data.recuperator.T_h-273.15;

% Expander
s_exp = cycle_data.expander.s/1000;
T_exp = cycle_data.expander.T-273.15;

% Pump
s_pump = cycle_data.pump_f.s/1000;
T_pump = cycle_data.pump_f.T-273.15;


%% Plot the T-s diagram
% Prepare the figure
fig = figure(); ax_fig = gca; 
hold on; box on;
pbaspect([1.15 1 1])
% axis square;

% Define the axis limits
TT = T_cycle(1:end);
ss = s_cycle(4:9);
T_minplot = min(0.90*T_min,min(TT)-(max(TT)-min(TT))/8);
T_maxplot = max(1.10*T_crit,max(TT)+(max(TT)-min(TT))/8);
if T_min < T_crit
    [~, s_sat] = sat_line(working_fluid,T_min+273.15,T_crit+273.15,'T','S',200); s_sat = s_sat/1000;
    s_minplot = min(s_sat(1)-(s_sat(end)-s_sat(1))/8,min(ss)-(max(ss)-min(ss))/8);
    s_maxplot = max(1.2*max(s_sat),max(ss)+(max(ss)-min(ss))/8);
else
    s_minplot = min(ss)-(max(ss)-min(ss))/8;
    s_maxplot = max(ss)+(max(ss)-min(ss))/8;
end
axis([s_minplot s_maxplot T_minplot T_maxplot])

% Label de axes
font_size = 12;
xlabel({' ';'$s$ -- Entropy (kJ/kg$\,$K)'},'FontSize',font_size);
ylabel({'$T$ -- Temperature ($^\circ$C)';' '},'FontSize',font_size);
ax_fig.YAxis.TickLabelFormat = '%.0f';
ax_fig.XAxis.TickLabelFormat = '%.2f';

% Plot theaturation line
[T_sat, s_sat] = sat_line(working_fluid,T_trip+273.15,T_crit+273.15,'T','S',200);
T_sat = T_sat-273.15; s_sat = s_sat/1000;
plot(s_sat,T_sat,'k','LineWidth',0.5)

% Plot the thermodynamic trajectories in the components
plot(s_pump,T_pump,'k')
plot(s_exp,T_exp,'k')
plot(s_rec_h,T_rec_h,'k')
plot(s_rec_c,T_rec_c,'k')
plot(s_cond_h,T_cond_h,'k')
plot(s_evap_c,T_evap_c,'k')
plot(s_cond_h,T_cond_h,'k')
plot(s_evap_c,T_evap_h,'r')
plot(s_cond_h,T_cond_c,'b')

% Plot the cycle points
plot(s_cycle(4:9),T_cycle(4:9),'ko','MarkerFaceColor','k','MarkerSize',2)
plot(s_cycle([7,6,6]),T_cycle([1,2,3]),'ro','MarkerFaceColor','r','MarkerSize',2)
plot(s_cycle([4,4,9]),T_cycle([10,11,12]),'bo','MarkerFaceColor','b','MarkerSize',2)
plot(s_cycle([6,6]),T_cycle([2,3]),'r-')
plot(s_cycle([4,4]),T_cycle([10,11]),'b-')

% Label the cycle points
p1 = '$\quad 1$';   text(s_cycle(7),T_cycle(1),     p1,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p2 = '$2 \quad$';   text(s_cycle(6),T_cycle(2)-3,   p2,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p3 = '$3 \quad$';   text(s_cycle(6),T_cycle(3)+3,   p3,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p4 = '$4 \quad$';   text(s_cycle(4),T_cycle(4)-3,   p4,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p5 = '$5 \quad$';   text(s_cycle(5),T_cycle(5)+3,   p5,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p6 = '$6 \quad$';   text(s_cycle(6),T_cycle(6)+6,   p6,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p7 = '$\quad 7$';   text(s_cycle(7),T_cycle(7),     p7,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p8 = '$\quad 8$';   text(s_cycle(8),T_cycle(8)+3,   p8,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p9 = '$\quad 9$';   text(s_cycle(9),T_cycle(9)-3,   p9,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p10 = '$10 \quad$'; text(s_cycle(4),T_cycle(10)-3,  p10, 'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p11 = '$11 \quad$'; text(s_cycle(4),T_cycle(11)+3,  p11, 'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p12 = '$\quad 12$'; text(s_cycle(9),T_cycle(12),    p12, 'HorizontalAlignment', 'left',  'FontSize', font_size-2)

% Save the figure
name = 'Ts_diagram_';
if save == 1
    saveas(fig,fullfile(my_path,[name,my_filename,'.pdf']),'pdf')
elseif save == 2
%     saveas(fig,fullfile(my_path,[name,my_filename]),'fig')
%     export_fig(fig,fullfile(my_path,[name,my_filename]),'-png','-r1000')
%     export_fig(fig,fullfile(my_path,[name,my_filename]),'-eps','-painters')
    export_fig(fig,fullfile(my_path,[name,my_filename]),'-pdf','-painters')
elseif save ~= 0
    error('Choose a valid saving option')
end

end
