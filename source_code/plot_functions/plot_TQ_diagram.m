function [] = plot_TQ_diagram(cycle_data,my_filename,my_path,save)
%% Load the trajectories in the components
% Evaporator
h_evap_h = cycle_data.evaporator.h_h/1000;
T_evap_h = cycle_data.evaporator.T_h-273.15;
T_evap_c = cycle_data.evaporator.T_c-273.15;

% Condenser
h_cond_h = cycle_data.condenser.h_h/1000;
T_cond_h = cycle_data.condenser.T_h-273.15;
T_cond_c = cycle_data.condenser.T_c-273.15;

% Recuperator
h_rec_h = cycle_data.recuperator.h_h/1000;
T_rec_h = cycle_data.recuperator.T_h-273.15;
T_rec_c = cycle_data.recuperator.T_c-273.15;

% Compute the heat flow rates
m_h = cycle_data.mass_flows.m_h;
m_f = cycle_data.mass_flows.m_f;
q_evap = (h_evap_h-h_evap_h(1))*m_h;
q_cond = (h_cond_h-h_cond_h(1))*m_f;
q_rec = (h_rec_h-h_rec_h(1))*m_f;


%% Plot the T-Q diagram
% Prepare the figure
fig = figure(); ax_fig = gca; 
hold on; box on; axis square;

% Label de axes
font_size = 12;
xlabel({' '; '$\dot{Q}$ -- Heat transfer (kW)'},'FontSize',font_size)
ylabel({'$T$ -- Temperature ($^\circ$C)';' '},'FontSize',font_size)
ax_fig.XAxis.TickLabelFormat = '%.0f';
ax_fig.YAxis.TickLabelFormat = '%.0f';

plot(0,0,'r')
plot(0,0,'k')
plot(0,0,'b')

plot(q_cond,T_cond_c,'b')
plot(q_cond,T_cond_h,'k')
plot((q_cond(end)+q_rec),T_rec_h,'k')
plot((q_cond(end)+q_rec),T_rec_c,'k')
plot((q_cond(end)+q_rec(end)+q_evap),T_evap_h,'r')
plot((q_cond(end)+q_rec(end)+q_evap),T_evap_c,'k')

x_lim = ax_fig.XLim;
y_lim = ax_fig.YLim;

plot([q_cond(end) q_cond(end)],[-10000 10000],'k','LineWidth',0.50)
plot([q_cond(end)+q_rec(end) q_cond(end)+q_rec(end)],[-10000 10000],'k','LineWidth',0.50)
plot([q_cond(end)+q_rec(end)+q_evap(end) q_cond(end)+q_rec(end)+q_evap(end)],[-10000 10000],'k','LineWidth',0.50)

ax_fig.XLim = x_lim;
ax_fig.YLim = y_lim;

% Save the figure
name = 'TQ_diagram';
if save == 1
    saveas(fig,fullfile(my_path,[my_filename,'_',name,'.pdf']),'pdf')
elseif save == 2
%     saveas(fig,fullfile(my_path,[my_filename,'_',name]),'fig')
%     export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-png','-r1000')
%     export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-eps','-painters')
    export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-pdf','-painters')
elseif save ~= 0
    error('Choose a valid saving option')
end

end
