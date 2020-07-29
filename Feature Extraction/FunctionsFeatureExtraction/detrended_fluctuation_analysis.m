function [n_win_tested, fn, alpha1, alpha2] = detrended_fluctuation_analysis(time, sig, plotFigure)

% DFA Detrended fluctuation analysis
% Method for determining the statistical self-affinity of a signal.
% Useful in revealing the extent of long-range correlations in time series
% https://www.physionet.org/physiotools/dfa/

% DFA is used to quantify the fractal scaling properties of heart rate
% signals of short interval. Describes the fractal-like correlation 
% properties of R–R interval data (Hoshi2013)

%   Calculates the DFA of a signal and it's scaling exponent alpha.
%   Input:
%       - time: time (or x values of signal)
%       - sig: signal data (or y values of signal)
%       - varargin: Pass in name-value pairs to configure advanced options:
%           - n_min: Minimal DFA block-size (default 4)
%           - n_max: Maximal DFA block-size (default 64)
%           - n_incr: Increment value for n (default 2)
%           - alpha1_range: Range of block size values to use for calculating the alpha1 scaling
%             exponent. Default: [4, 15].
%           - alpha2_range: Range of block size values to use for calculating the alpha2 scaling
%             exponent. Default: [16, 64].
%   Output:
%       - n_win_tested: block sizes (x-axis of DFA)
%       - fn: DFA value for each block size n
%       - alpha: Exponential scaling factor


% DFA plot is not strictly linear but rather consisted of two distinct
% regions of different slopes separated at a break point suggesting that 
% there is a short range scaling exponent (alpha1) over periods of 4 to 11 
% beats (or 4 to 13), and a long-range exponent (alpah2), over longer periods 
% (larger than 11 beats) (Hoshi2013)

n_min = 4;
n_max = 64;
n_incr = 2;
alpha1_range = [4, 11]; % from Hoshi2013 and Costa2017
alpha2_range = [11, 64]; % from Hoshi2013 and Costa2017

alpha1 = NaN;
alpha2 = NaN;


%% DFA

% Integrate the signal without mean
% Given a bounded time series x_t of length N, integration or summation
% first converts this into an unbounded process X_t:
nni_int = cumsum(sig - mean(sig));
% nni_int is called cumulative sum or profile

N = length(nni_int);

fn = NaN(n_max, 1);
n_win_tested = n_min:n_incr:n_max;

% The n-th order polynomial regressor in the DFA family is typically 
% denoted as DFAn, with unlabeled DFA often referring to DFA1 (Bryce2012).

for nn = n_win_tested % process repeated over a range of different window sizes nn 
    % to characterize the relationship between fn, the average fluctuation,
    % and the window size, n_win_tested
    % Typically, fn will increase with window size. 
    
    % Calculate the number of windows we need for the current n
    num_win = floor(N/nn);

    % Break the signal into num_win windows of n samples each
    sig_windows = reshape(nni_int(1:nn*num_win), nn, num_win);
    t_windows  = reshape(time(1:nn*num_win), nn, num_win);
    sig_regressed = zeros(size(sig_windows));
%     sig_regressed2 = zeros(size(sig_windows));
    
    % Perform linear regression in each window
    % least squares line is fit to the data 
    % (representing the trend in that box)
    for ii = 1:num_win
        y = sig_windows(:, ii);
        x = [ones(nn, 1), t_windows(:, ii)];
        slope = x\y;
        yn = x * slope;
        sig_regressed(:, ii) = yn;
%         fitted = polyfit(t_windows(:, ii)', sig_windows(:, ii)', 1);
%         sig_regressed2(:, ii) = t_windows(:, ii).*fitted(1)+fitted(2);
    end

    % Calculate the fluctuation F(n), the value of the DFA for the current n
    % in other words the root-mean-square deviation from the trend
    fn(nn) = sqrt ( 1/N * sum((sig_windows(:) - sig_regressed(:)).^2) );
    % 1st the integrated time-series (nni_int) is detrended, by subtracting
    % the local trend yn for each n_win_tested
end


% Find the indices of all the DFA values we calculated
fn = fn(n_win_tested);
n_win_tested = n_win_tested';

if fn==0
    disp('problem')
end


fn = fn(fn~=0);
n_win_tested = n_win_tested(fn~=0);
% last_n_win_tested = n_win_tested(end);

%% Scaling exponent, alpha

% Fit a linear regression model to the log-log DFA in each alpha range
% to quantify the fractal scaling properties of RR intervals.
fn_log = log10(fn);
n_log = log10(n_win_tested);

% Find DFA values in each of the alpha ranges
alpha1_idx = find(n_win_tested >= alpha1_range(1) & n_win_tested <= alpha1_range(2));
alpha2_idx = find(n_win_tested >= alpha2_range(1) & n_win_tested <= alpha2_range(2));

alpha_idx_vec = {alpha1_idx; alpha2_idx};

degree = 1;
fit_vec = zeros(2,2);
for ii = 1:2
    
    alpha_idx = cell2mat(alpha_idx_vec(ii));
    
    if numel(alpha_idx)<2
       continue 
    end
    
    % Least-squares fit polynomial coefficients
    fit_vec(ii,:) = polyfit(n_log(alpha_idx), fn_log(alpha_idx), degree); 
    
    % Save the slopes of the lines = scaling exponents
    alpha_val = fit_vec(ii,1);
    % alpha is a measure of correlation in the noise and is simply an
    % estimate of the Husrt exponent H (Bryce2012)

    % Interpretation of alpha values
    eval(['alpha' num2str(ii) ' = alpha_val;'])
%     disp(['alpha' num2str(ii) ' = ' num2str(alpha_val)])
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
    

    figure()
    lw = 3.8; ls = ':';
    
    loglog(n_win_tested, fn, 'ko', 'MarkerSize', 7); % h1 = 
    hold on
    grid on
    axis tight
    
    alpha_legend = {};
    color_vec = {'b','r'};
    for ii = 1:2
        
        alpha_idx = cell2mat(alpha_idx_vec(ii));
        
        if numel(alpha_idx)<2
            continue
        end
        
        % Plot alpha1 line
        alpha_line = 10.^(fit_vec(ii,1) * n_log(alpha_idx) + fit_vec(ii,2));
        loglog(n_win_tested(alpha_idx), alpha_line, 'Color', color_vec{ii}, 'LineStyle', ls, 'LineWidth', lw);
        
        % Calculate fit R-square, to include in legend
        detailed_legend = true;
        if detailed_legend
            C = corrcoef(fn(alpha_idx), alpha_line);
            R2 = C(1,2)^2;
            
            alpha_legend{ii} = ['$\alpha_' num2str(ii) ' = ' ...
                num2str(round(fit_vec(ii,1)*1000)/1000)...
                ' \, (R^2 = ' num2str(round(R2*1000)/1000) ')$'];
        else
            alpha_legend{ii} = ['$\alpha_' num2str(ii) ' = '...
                num2str(round(fit_vec(ii,1)*1000)/1000) '$'];
            
        end
        
    end
    

    xlabel('$log(n)$'); ylabel('$log\big[ F(n) \big]$');
    set(gca,'XTick', 2.^(1:15)); % Set ticks at powers of two

    legend(['DFA', alpha_legend], 'Location', 'northwest');
%     uistack(h1, 'top');

    
end


end