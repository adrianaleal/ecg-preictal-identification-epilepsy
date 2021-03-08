function [] = copy_figure(old_fig_handle, new_plot_handle)

% Identify axes to be copied:
axes_to_be_copied = findobj(old_fig_handle,'type','axes');


% Identify the children of this axes:
children_to_be_copied = get(axes_to_be_copied,'children');

copyobj(children_to_be_copied,new_plot_handle);

xlabel(old_fig_handle.CurrentAxes.XLabel.String)
ylabel(old_fig_handle.CurrentAxes.YLabel.String)
title(old_fig_handle.CurrentAxes.Title.String)
hLines = findobj(new_plot_handle, 'Type', 'Line');
legend([hLines(end) hLines(1)], old_fig_handle.CurrentAxes.Legend.String, ...
    'Location', 'Best')
set(new_plot_handle, 'YScale', old_fig_handle.CurrentAxes.YScale)
set(new_plot_handle, 'XScale', old_fig_handle.CurrentAxes.XScale)
axis tight
grid on
box on
% Close the old figure 
close(old_fig_handle); 


end