function [SD1, SD2, SD1toSD2] = poincare_plot(RRI_signal_segment, plotFigure)

% Read Shaffer2017 for interpretation of the results

% this code was adapted from function poincare.m in 
% https://github.com/physiozoo/mhrv

% Inputs:
% - RRI_signal_segment (double): RR interval series
% - plotFigure (double): flag to plot figure

% Outputs:
% - SD1: Poincare plot SD1 descriptor (std. dev. of intervals along
%        the line perpendicular to the line of identity).
% - SD2: Poincare plot SD2 descriptor (std. dev. of intervals along
%        the line of identity).
% - SD1toSD2: ratio of SD1 to SD2



RRi = RRI_signal_segment(1:end-1); % RR(i)
RRii = RRI_signal_segment(2:end); % RR(i+1)


% from the geometrical analysis (in Brennan2001):
% x1 = cos(pi/4)*RR(i) - sin(pi/4)*RR(i+1)
% <=> x1 = 1/sqrt(2)*RR(i) - 1/sqrt(2)*RR(i+1)

% x2 = sin(pi/4)*RR(i) + cos(pi/4)*RR(i+1)
% <=> x2 = 1/sqrt(2)*RR(i) + 1/sqrt(2)*RR(i+1)




% SD1:
SD1 = std(RRi-RRii)/sqrt(2);

% SD2:
SD2 = std(RRi+RRii)/sqrt(2);

% SDRR
% SDRR = sqrt(SD1^2+SD2^2)/sqrt(2);

SD1toSD2 = SD1/SD2;

if plotFigure
    
    
    %% Rotate input
    alpha = -pi/4;

    % Rotate the data to the new coordinate system
    % _new is the square-filtered data in the new coordinate system
    RR_interval_signal_segment_rotated = rotation_matrix(alpha) * [RRi; RRii];
    RRi_rotated = RR_interval_signal_segment_rotated(1,:);
    RRii_rotated = RR_interval_signal_segment_rotated(2,:);

    % Calculate standard deviation along the new axes
%     sd1 = sqrt(var(RRii_rotated));
%     sd2 = sqrt(var(RRi_rotated));

    sd1_factor = 1;
%     sd2_factor = sd1_factor;
    %% Fit ellipse
    % For fitting the ellipse we're using the _new vectors because we don't wan't any non-physiological
    % intervals to influence the ellise and SD1/2 metrics.
    
    % Ellipse radii
    r_x = sd1_factor * SD2;
    r_y = sd1_factor * SD1;
    
    % Ellipse center
    c_x = mean(RRi_rotated);
    c_y = mean(RRii_rotated);
    
    % Ellipse parametric equation
    t = linspace(0, 2*pi, 200);
    xt = r_x * cos(t) + c_x;
    yt = r_y * sin(t) + c_y;
    
    % Rotate the ellipse back to the old coordinate system
    ellipse_original_axis = rotation_matrix(-alpha) * [xt; yt];
    
    %% Lines for ellipse axes
    ellipse_center_rotated = [c_x; c_y];
    
    % Create the lines in the new coordinate system
    sd1_line_new = [0, 0; -SD1, SD1] + [ellipse_center_rotated, ellipse_center_rotated];
    sd2_line_new = [-SD2, SD2; 0, 0] + [ellipse_center_rotated, ellipse_center_rotated];
    
    % Rotate back to the old system
    sd1_line_old = rotation_matrix(-alpha) * sd1_line_new;
    sd2_line_old = rotation_matrix(-alpha) * sd2_line_new;
    
    %% Plotting
    figure();
    plot(RRi,RRii,'*'), hold on
    plot(RRi-RRii, 'or')
    plot(RRi+RRii, 'og')
    plot(ellipse_original_axis(1,:), ellipse_original_axis(2,:),'k--', 'LineWidth',1.5);
    h_sd1 = plot(sd1_line_old(1,:), sd1_line_old(2,:), 'r-', 'LineWidth', 3);
    h_sd2 = plot(sd2_line_old(1,:), sd2_line_old(2,:), 'g-', 'LineWidth', 3);
    grid on
    axis tight
    % [xmin xmax ymin ymax] 
    val = 30;
    axis([min(sd2_line_old(1,:))-val max(sd2_line_old(1,:))+val ...
        min(sd2_line_old(2,:))-val max(sd2_line_old(2,:))+val])
    
    % if latex:
    % xlabel('$RR_i$ (ms)')
    % ylabel('$RR_{i+1}$ (ms)')
    % legend([h_sd1,h_sd2], {['$SD_1 = ' num2str(round(SD1*1e4)/1e4) ...
    %     '$ (short-term HRV)'], ['$SD_2 = ', num2str(round(SD2*1e4)/1e4) ...
    %     '$ (long-term HRV)']});
    
    xlabel('RR_i (ms)')
    ylabel('RR_{i+1} (ms)')
    legend([h_sd1,h_sd2], {['SD_1 = ' num2str(round(SD1*1e4)/1e4) ...
        ' (short-term HRV)'], ['SD_2 = ', num2str(round(SD2*1e4)/1e4) ...
        ' (long-term HRV)']});
    
    title('RR Interval Poincare Plot')

end

end

function rotation_mat = rotation_matrix(theta)
%% Helper functions

% Creates a 2D rotation matrix

rotation_mat = [cos(theta), -sin(theta); sin(theta), cos(theta)];
end
