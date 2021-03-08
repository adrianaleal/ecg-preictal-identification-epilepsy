function [rqa_stat, distance_mat, attractor] = recurrenceAnalysis(input_signal, ...
    tau, eDim, plotFigure, attractor)

% Input:
%   - input_signal: signal data
%   - tau: time delay
%   - eDim: embedding dimension
%   - attractor: reconstructed phase space
%   - plotFigure: flag to plot figure

% Output:
%   - rqa_stat: recurrence quantification analysis statistics (REC DET
%     LMAX L ENT LAM TT)
%   - distance_mat: square matrix containing the pairwise Euclidean 
%     distance between all samples of the attractor
%   - attractor: reconstructed phase space

% Sources: 
% - Marwan2007: MARWAN, N., CARMENROMANO, M., THIEL, M. & KURTHS, J. 
%   Recurrence plots for the analysis of complex systems. Phys. Reports 
%   438, 237–329, DOI: 10.1016/j.physrep.2006.11.001 (2007).
% - Marwan2015: Marwan, N. & Webber, C. L. Mathematical and Computational 
%   Foundations of Recurrence Quantifications. In Webber, J. C. & Marwan, 
%   N. (eds.) Recurrence Quantification Analysis. Understanding Complex
%   Systems., chap. 1, 3–43, DOI: 10.1007/978-3-319-07155-8 1 (Springer, 
%   Cham, 2015).
% - Billeci, L., Marino, D., Insana, L., Vatti, G. & Varanini, M. 
%   Patient-specific seizure prediction based on heart rate variability and 
%   recurrence quantification analysis. PLOS ONE 13, 1–21, DOI: 
%   10.1371/journal.pone.0204339 (2018).


if nargin<4
    plotFigure = 0;
end

if nargin<5
    attractor = getAttractor(input_signal, tau, eDim, plotFigure);
end


%% get color recurrence plot
% square matrix of [N x N] containing the pairwise Euclidean distance
% between all samples of the attractor, N, in the eDim-dimensional space:
distance_mat = pdist2(attractor,attractor);

if plotFigure==1
    figure()
    subplot(121)
    imagesc(distance_mat)
    set(gca, 'YDir', 'normal')
    colormap Jet
    h = colorbar;
    ylabel(h, {'Distance between subsequent', 'points in the attractor'})
    % axis image
    xlabel('RRI series index')
    ylabel('RRI series index')
    title('Distance matrix')
end


%% get the recurrence plot

% get the maximal phase space diameter = |xMax – xMin|:
maxDiam = abs(max(attractor(:))-min(attractor(:)));

% get epsilon = 0.1|xMax – xMin|:
epsilon = 0.1 * maxDiam;

% get the binary (black-white) recurrence matrix after defining threshold 
% epsilon:
bin_recurr_mat = distance_mat<=epsilon;
[row_ind_recurrence_points, col_ind_recurrence_points] = find(bin_recurr_mat);
ind_recurrence_points = [row_ind_recurrence_points, col_ind_recurrence_points];

if plotFigure==1
    subplot(122)
    plot(row_ind_recurrence_points, col_ind_recurrence_points,'k.', ...
        'MarkerSize',2)
    xlim([0, size(distance_mat,1)])
    ylim([0, size(distance_mat,2)])
    % axis image
    xlabel('RRI series index')
    ylabel('RRI series index')
    title('Recurrence Plot')
end


% Recurrence quantification analysis
% rqa_stat - RQA statistics - [REC DET LMAX L ENT LAM TT]
% The main advantage of the recurrence quantification analysis is that it 
% can provide useful information even for short and non-stationary data

rqa_stat = recurrence_quantification_analysis(bin_recurr_mat, ...
    ind_recurrence_points);


end