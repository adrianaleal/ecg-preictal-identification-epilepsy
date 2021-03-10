function [segment_size_vec, n_segment_vec, F, alpha_vec] = ...
    monofractal_detrended_fluctuation_analysis(input_signal, scale_range, ...
    plotFigure)

% DFA Detrended fluctuation analysis
% Method for determining the statistical self-affinity of a signal.
% Useful in revealing the extent of long-range correlations in time series
 
% This function computes the DFA of a signal and it's scaling exponent 
% alpha.

% Input:
% - time: time (or x values of signal)
% - input_signal: signal data (or y values of signal)
% - scale_range: Range of segment size values to use for
%                calculating the alpha scaling exponent. The input should 
%                be a matrix of the number of alpha ranges (rows) by two 
%                columns corresponding to the minimum and maximum number of 
%                samples in a segment, respectively.

% Output:
% - segment_size_vec: vector containing the different segments 
%                     (x-axis of DFA)
% - n_segment_vec: vector containing the number of segments for each
%                  segment size
% - fn: DFA value for each segment size
% - alpha_vec: Exponential scaling factors corresponding to the
%              ranges defined in alpha_range matrix


% READ:
% - https://www.physionet.org/physiotools/dfa/
% - Ihlen2012: Ihlen, E. A. F. Introduction to Multifractal Detrended 
%   Fluctuation Analysis in Matlab. Front. Physiol. 3, 141, DOI:
%   10.3389/fphys.2012.00141 (2012).
% - Bryce2012: Bryce, R. M. & Sprague, K. B. Revisiting detrended 
%   fluctuation analysis. Sci. Reports 2, 315, DOI: 10.1038/srep00315 
%   (2012).
% - MoralesMartinez2021: Morales Martínez, J. L., Segovia-Domínguez, I., 
%   Rodríguez, I. Q., Horta-Rangel, F. A. & Sosa-Gómez, G. A modified
%   Multifractal Detrended Fluctuation Analysis (MFDFA) approach for 
%   multifractal analysis of precipitation. Phys. A: Stat. Mech. its Appl. 
%   565, 125611, DOI: 10.1016/j.physa.2020.125611 (2021).
% - Hoshi2013: Hoshi, R. A., Pastre, C. M., Vanderlei, L. C. M. & Godoy, M.
%   F. Poincaré plot indexes of heart rate variability: Relationships with 
%   other nonlinear variables. Auton. Neurosci. 177, 271–274, DOI: 
%   10.1016/j.autneu.2013.05.004 (2013).
% - Costa2017: Costa, M. D., Davis, R. B. & Goldberger, A. L. Heart Rate 
%   Fragmentation: A New Approach to the Analysis of Cardiac Interbeat 
%   Interval Dynamics. Front. Physiol. 8, 1–13, DOI: 
%   10.3389/fphys.2017.00255 (2017).
% - Voss2009: Voss, A., Schulz, S., Schroeder, R., Baumert, M. & Caminal, 
%   P. Methods derived from nonlinear dynamics for analysing heart rate 
%   variability. Philos. Transactions Royal Soc. A: Math. Phys. Eng. Sci. 
%   367, 277–296, DOI: 10.1098/rsta.2008.0232 (2009).

% Limitations: at least 8000 data points should be used; monofractal method;
% (Bryce2012, Ihlen2012, MoralesMartinez2021, Voss2009)

% Limitations on ECG: normal-to-normal interbeat intervals are required; 
% dependency on editing ectopic beats. Voss2009



%% CONSIDERATIONS ON ECG DATA:
% DFA is used to quantify the fractal scaling properties of heart rate
% signals of short interval. Describes the fractal-like correlation 
% properties of R–R interval data (Hoshi2013)

% DFA plot is not strictly linear but rather consisted of two distinct
% regions of different slopes separated at a break point suggesting that 
% there is a short range scaling exponent (alpha1) over periods of 4 to 11 
% beats (or 4 to 13), and a long-range exponent (alpha2), over longer 
% periods (larger than 11 beats) (Hoshi2013)

% scale_range = [4 11; 11 64]; % from Hoshi2013 and Costa2017

%%

N = length(input_signal);

if nargin == 1
    scale_range = [4 11; 11 64];
    plotFigure = 0;
end

if nargin == 2
    plotFigure = 0;
end


% get the minimum number of samples in a segment
min_subsequence_length = min(scale_range(:));
% Ihlen2012: "The minimum segment size larger than 10 samples is a "rule of
% thumb" for the computation of RMS."


% get the maximum number of samples in a segment
max_subsequence_length = max(scale_range(:));
% Ihlen2012: "A maximum segment size below 1/10 of the sample size of the 
% time series will provide at least 10 segments in the computation of F."


% define the number of divisions:
n_scales = 50;
exponents = linspace(log2(min_subsequence_length),log2(max_subsequence_length),n_scales); 
segment_size_vec = unique(floor(2.^exponents));
n_scales = length(segment_size_vec); % update the number of divisions
% Ihlen2012: "It is favorable to have a equal spacing between scales when 
% they are represented in the log-log DFA."



%% DFA


% Convert a noise like time series into a random walk like time series:
rw_input_signal = cumsum(input_signal - mean(input_signal));
% rw_input_signal = input_signal;
% OR integrate the signal without mean
% OR given a bounded time series of length N, integration or summation
% first converts this into an unbounded process 
% OR rw_input_signal is called cumulative sum or profile


F = NaN(n_scales, 1);

n_segment_vec = NaN(n_scales,1);

% define the order of the polynomial fit
m = 1; % if the polinomial trend is linear
% m = 2; % if the polinomial trend is quadratic
% m = 3; % if the polinomial trend is cubic


% The n-th order polynomial regressor in the DFA family is typically 
% denoted as DFAn, with unlabeled DFA often referring to DFA1 (Bryce2012).

indexes = 1:N;

for n = 1:n_scales
    % process repeated over a range of different segments lengths
    % to characterize the relationship between fn, the average fluctuation,
    % and the segment size
    % Typically, Fn will increase with segment size. 
    
    
    segment_size = segment_size_vec(n); % size of the segments
    n_segment = floor(N/segment_size); % number of segments
    n_segment_vec(n) = n_segment;

    RMS = zeros(n_segment,1);

    
    % Perform local detrending (or linear regression) in each segment
    % least squares line is fit to the data 
    % (representing the trend in that segment)
    for jj = 1:n_segment
        % Break the signal into n_segment segments of
        % segment_size each
    
        ind_start = ((jj-1)*segment_size)+1; 
        ind_stop = jj*segment_size; 
        indexes_segment = (ind_start:ind_stop)';
        signal_segment = rw_input_signal(indexes_segment)';
        
        % (1) Fit a polynomial trend to each segment        
        if m==1
            % OPTION1: lowest computational effort for a linear fit
            x = [ones(segment_size, 1), indexes_segment];
            slope = x\signal_segment; % x/y
            yn = x * slope;
            segments_trend = yn;
        else
            % OPTION2: for fitting polynomials of higher order
            % get the polinomial coefficients C used to create the
            % polynomial trend segments_trend(:, jj):
            C = polyfit(indexes_segment, signal_segment, m);
            segments_trend = polyval(C, indexes_segment);
            % OR:
            % segments_trend(:, jj) = indexes_segments(:, jj).*C(1)+C(2);
        end
               
        % (2) Compute the local fluctuation or RMS around a trend 
        % segments_trend(:, jj), for each segment
        residual_variation = signal_segment-segments_trend;
        RMS(jj) = sqrt(mean(residual_variation.^2));
        
    end
    
    % Compute the overall fluctuation or overall RMS OR scaling function:
    F(n) = sqrt(mean(RMS.^2));
    % Calculate the fluctuation F, the value of the DFA for the current
    % segment size
    % in other words the root-mean-square deviation from the trend

end



if F==0
    disp('problem')
end


F = F(F~=0);
segment_size_vec = segment_size_vec(F~=0);
segment_size_vec = segment_size_vec';

% last_n_win_tested = n_win_tested(end);

%% Scaling exponent, alpha

% Fit a linear regression model to the log-log DFA in each alpha range
% to quantify the fractal scaling properties of the input input_signal.

% Ihlen2012: DFA indentifies the monofractal structure of a time series as
% the power law relation between the overall RMS computed for multiple 
% scales (defined in each alpha range)


% get the number of scale ranges:
n_scale_range = size(scale_range,1);

% initialize output result
alpha_vec = NaN(n_scale_range,1);


m = 1; % define the order of the polynomial fit
fit_vec = zeros(n_scale_range,2);

% Find DFA values in each of the alpha ranges
alpha_idx_vec = cell(n_scale_range,1);
for aa = 1:n_scale_range
    alpha_idx_vec{aa} = find(segment_size_vec >= scale_range(aa,1) & segment_size_vec <= scale_range(aa,2));
    
    alpha_idx = cell2mat(alpha_idx_vec(aa));
    
    if numel(alpha_idx)<2
       continue 
    end
    
    % Least-squares fit polynomial coefficients
    fit_vec(aa,:) = polyfit(log10(segment_size_vec(alpha_idx)), ...
        log10(F(alpha_idx)), m); 
    
    % Save the slopes of the lines = scaling exponents
    alpha_val = fit_vec(aa,1);
    alpha_vec(aa) = alpha_val;
    % alpha is a measure of correlation in the noise and is simply an
    % estimate of the Hurst exponent H (Bryce2012)

    % Interpretation of alpha values
    eval(['alpha' num2str(aa) ' = alpha_val;'])
%     disp(['alpha' num2str(aa) ' = ' num2str(alpha_val)])
%     if (0 < alpha_val) && (alpha_val < 0.5)
%         disp('Large and small values of the time series are more likely to alternate!')
%     elseif alpha_val==0.5
%         disp('White noise')
%     elseif (0.5 < alpha_val) && (alpha_val < 1.0)
%         disp('Persistent long-range power-law correlations such that a large')
%         disp('(compared to the average) interbeat interval is more likely to be')
%         disp('followed by large interval and vice versa')
%     elseif alpha_val == 1
%         disp('1/f noise, pink noise')
%     elseif alpha_val > 1
%         disp('Correlations exist but cease to be of a power-law form')
%         if alpha_val==1.5
%             disp('Brownian noise or random walk')
%         end
%     end

end

%% Plot
% log-log graph of F(n) against n

if plotFigure
    
    lw = 3.8; ls = ':';
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(221)
    plot(1:N, input_signal, 'b'), hold on
    plot(1:N, rw_input_signal, 'r'), hold off
    axis tight
    ylabel('Amplitude (measurement units)','FontSize',14);
    xlabel('time (sample number)','FontSize',14);
    legend('Noise like input signal', 'Random walk like input signal','FontSize',14)
    
    subplot(222)
    % plot of the local trends for a given segment size
    plot(1:N, rw_input_signal, 'r'), hold on
    
    segment_size_plot = 64;
    index = find(segment_size_vec>=segment_size_plot, 1);
    segment_size_plot = segment_size_vec(index);
    n_segment_plot = n_segment_vec(index);
    
    indexes_segments = reshape(indexes(1:n_segment_plot*segment_size_plot), ...
        segment_size_plot, n_segment_plot);
    
        
    RMS = zeros(segment_size_plot, n_segment_plot);
    for jj = 1:n_segment_plot
        
        signal_segment = rw_input_signal(indexes_segments(:, jj))';
        
        C = polyfit(indexes_segments(:, jj), signal_segment, m);
        segments_trend = polyval(C,indexes_segments(:, jj));  
        
        residual_variation = signal_segment-segments_trend;
        RMS(:,jj) = sqrt(mean(residual_variation.^2)); 
        
        mat = [segments_trend, segments_trend+RMS(1,jj), ...
            segments_trend-RMS(1,jj)];
        
        pplot = plot(indexes_segments(:,jj), mat, 'b');
        set(pplot(1),'LineStyle','--');
        set(pplot(2),'LineStyle','-'); 
        set(pplot(3),'LineStyle','-');
        
    end
    
    F_plot = sqrt(mean(RMS(1,:).^2));
    
    x = indexes_segments(end, :);
    y = [min(rw_input_signal) max(rw_input_signal)];
    line(repmat(x,2,1),y, 'Color','k','LineStyle','--')
    hold off
    axis tight
    ylabel('Amplitude (measurement units)','FontSize',14);
    xlabel('time (sample number)','FontSize',14);
    legend('Random walk like input signal', ...
        'Local trend for each segment', '+/- local RMS', 'FontSize',14)
    title(['segment size = ' num2str(segment_size_plot) ' samples'])
    
    subplot(223)
    % plot of the local fluctuations for a given segment size
    plot(indexes_segments(:), RMS(:)), hold on
    plot(indexes_segments(:), repmat(F_plot,1, segment_size_plot*n_segment_plot)), hold off
    axis tight
    ylabel('Amplitude (measurement units)','FontSize',14);
    xlabel('time (sample number)','FontSize',14);
    legend('Local fluctuations, RMS, for each segment', ...
        'Overall fluctuation, F ', 'FontSize',14)
    title(['segment size = ' num2str(segment_size_plot) ' samples'])
    
    subplot(224)
    loglog(segment_size_vec, F, 'ko', 'MarkerSize', 7);
    hold on
    grid on
    axis tight
    
    alpha_legend = {};
    color_vec = {'b','r','m'};
    for aa = 1:n_scale_range
        
        alpha_idx = cell2mat(alpha_idx_vec(aa));
        
        if numel(alpha_idx)<2
            continue
        end
        
        % Plot alpha1 line
        alpha_reg_line = polyval(fit_vec(aa,:),log10(segment_size_vec(alpha_idx)));
        loglog(segment_size_vec(alpha_idx), 10.^(alpha_reg_line), ...
            'Color', color_vec{aa}, 'LineStyle', ls, 'LineWidth', lw);
        
        
        % Calculate fit R-square, to include in legend
        detailed_legend = true;
        if detailed_legend
            C = corrcoef(F(alpha_idx), alpha_reg_line);
            R2 = C(1,2)^2;
            
            % if latex
            % alpha_legend{aa} = ['$\alpha_' num2str(aa) ' = ' ...
            %     num2str(round(fit_vec(aa,1)*1000)/1000)...
            %     ' \, (R^2 = ' num2str(round(R2*1000)/1000) ')$'];
            
            if n_scale_range==1
                alpha_legend{aa} = ['\alpha = ' ...
                    num2str(round(fit_vec(aa,1)*1000)/1000)...
                    '  (R^2 = ' num2str(round(R2*1000)/1000) ')'];
            else
                alpha_legend{aa} = ['\alpha_' num2str(aa) ' = ' ...
                    num2str(round(fit_vec(aa,1)*1000)/1000)...
                    '  (R^2 = ' num2str(round(R2*1000)/1000) ')'];
            end
        else
            % if latex
            % alpha_legend{aa} = ['$\alpha_' num2str(aa) ' = '...
            %     num2str(round(fit_vec(aa,1)*1000)/1000) '$'];
            
            if n_scale_range==1
                alpha_legend{aa} = ['\alpha = '...
                    num2str(round(fit_vec(aa,1)*1000)/1000)];
            else
                alpha_legend{aa} = ['\alpha_' num2str(aa) ' = '...
                    num2str(round(fit_vec(aa,1)*1000)/1000)];
            end
        end
        
    end
    
    % xlabel('$log(n)$'); ylabel('$log\big[ F(n) \big]$'); % if latex:
    xlabel('log10(n)'); ylabel('log10[F(n)]');
    set(gca,'XTick', 2.^(1:15)); % Set ticks at powers of two

    legend(['DFA', alpha_legend], 'Location', 'northwest');
    
end


end