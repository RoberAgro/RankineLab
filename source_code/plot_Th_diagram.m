function [] = plot_Th_diagram(cycle_data,my_filename,my_path,save)
%% Load the cycle points
h_cycle = [cycle_data.cycle_states.h]/1000;
T_cycle = [cycle_data.cycle_states.T]-273.15;
working_fluid = cycle_data.fluids.working_fluid;
T_min = cycle_data.properties.T_min-273.15;
T_trip = cycle_data.properties.T_trip-273.15;
T_crit = cycle_data.properties.T_crit-273.15;


%% Load the trajectories in the components
% Evaporator
h_evap_c = cycle_data.evaporator.h_c/1000;
T_evap_c = cycle_data.evaporator.T_c-273.15;
T_evap_h = cycle_data.evaporator.T_h-273.15;

% Condenser
h_cond_h = cycle_data.condenser.h_h/1000;
T_cond_h = cycle_data.condenser.T_h-273.15;
T_cond_c = cycle_data.condenser.T_c-273.15;

% Recuperator
h_rec_c = cycle_data.recuperator.h_c/1000;
h_rec_h = cycle_data.recuperator.h_h/1000;
T_rec_c = cycle_data.recuperator.T_c-273.15;
T_rec_h = cycle_data.recuperator.T_h-273.15;

% Expander
h_exp = cycle_data.expander.h/1000;
T_exp = cycle_data.expander.T-273.15;

% Pump
h_pump = cycle_data.pump_f.h/1000;
T_pump = cycle_data.pump_f.T-273.15;


%% Plot the T-h diagram
% Prepare the figure
fig = figure(); ax_fig = gca; 
hold on; box on;
pbaspect([1.15 1 1])
% axis square;

% Define the axis limits
TT = T_cycle(1:end);
hh = h_cycle(4:9);
T_minplot = min(0.90*T_min,min(TT)-(max(TT)-min(TT))/8);
T_maxplot = max(1.10*T_crit,max(TT)+(max(TT)-min(TT))/8);
if T_min < T_crit
    [~, h_sat] = sat_line(working_fluid,T_min+273.15,T_crit+273.15,'T','H',200); h_sat = h_sat/1000;
    h_minplot = min(h_sat(1)-(h_sat(end)-h_sat(1))/8,min(hh)-(max(hh)-min(hh))/8);
    h_maxplot = max(1.2*max(h_sat),max(hh)+(max(hh)-min(hh))/8);
else
    h_minplot = min(hh)-(max(hh)-min(hh))/8;
    h_maxplot = max(hh)+(max(hh)-min(hh))/8;
end

axis([h_minplot h_maxplot T_minplot T_maxplot])

% Label de axes
font_size = 12;
xlabel({' ';'$h$ -- Enthalpy (kJ/kg)'},'FontSize',font_size);
ylabel({'$T$ -- Temperature ($^\circ$C)';' '},'FontSize',font_size);
ax_fig.YAxis.TickLabelFormat = '%.0f';
ax_fig.XAxis.TickLabelFormat = '%.0f';

% Plot theaturation line
[T_sat, h_sat] = sat_line(working_fluid,T_trip+273.15,T_crit+273.15,'T','H',200);
T_sat = T_sat-273.15; h_sat = h_sat/1000;
plot(h_sat,T_sat,'k','LineWidth',0.5)

% Plot the thermodynamic trajectories in the components
plot(h_pump,T_pump,'k')
plot(h_exp,T_exp,'k')
plot(h_rec_h,T_rec_h,'k')
plot(h_rec_c,T_rec_c,'k')
plot(h_evap_c,T_evap_c,'k')
plot(h_cond_h,T_cond_h,'k')
plot(h_evap_c,T_evap_h,'r')
plot(h_cond_h,T_cond_c,'b')

% Plot the cycle points
plot(h_cycle(4:9),T_cycle(4:9),'ko','MarkerFaceColor','k','MarkerSize',2)
plot(h_cycle([7,6,6]),T_cycle([1,2,3]),'ro','MarkerFaceColor','r','MarkerSize',2)
plot(h_cycle([4,4,9]),T_cycle([10,11,12]),'bo','MarkerFaceColor','b','MarkerSize',2)
plot(h_cycle([6,6]),T_cycle([2,3]),'r-')
plot(h_cycle([4,4]),T_cycle([10,11]),'b-')

% Label the cycle points
p1 = '$\quad 1$';   text(h_cycle(7),T_cycle(1),     p1,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p2 = '$2 \quad$';   text(h_cycle(6),T_cycle(2)-3,   p2,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p3 = '$3 \quad$';   text(h_cycle(6),T_cycle(3)+3,   p3,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p4 = '$4 \quad$';   text(h_cycle(4),T_cycle(4)-3,   p4,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p5 = '$5 \quad$';   text(h_cycle(5),T_cycle(5)+3,   p5,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p6 = '$6 \quad$';   text(h_cycle(6),T_cycle(6)+6,   p6,  'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p7 = '$\quad 7$';   text(h_cycle(7),T_cycle(7),     p7,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p8 = '$\quad 8$';   text(h_cycle(8),T_cycle(8)+3,   p8,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p9 = '$\quad 9$';   text(h_cycle(9),T_cycle(9)-3,   p9,  'HorizontalAlignment', 'left',  'FontSize', font_size-2)
p10 = '$10 \quad$'; text(h_cycle(4),T_cycle(10)-3,  p10, 'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p11 = '$11 \quad$'; text(h_cycle(4),T_cycle(11)+3,  p11, 'HorizontalAlignment', 'right', 'FontSize', font_size-2)
p12 = '$\quad 12$'; text(h_cycle(9),T_cycle(12),    p12, 'HorizontalAlignment', 'left',  'FontSize', font_size-2)

% Save the figure
name = 'Th_diagram_';
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
