function rqa_stat = recurrence_quantification_analysis(recurrence_plot, ...
    ind_recurrence_points, lmin)

% Input:
%   - recurrence_plot: binary recurrence matrix
%   - ind_recurrence_points: time delay
%   - lmin: minimum number of points that form vertical and diagonal lines

% Output:
%   - rqa_stat: recurrence quantification analysis statistics (REC DET
%     LMAX L ENT LAM TT)

if nargin<3 || isempty(lmin)
    lmin = 2;
end

% LOI: line of identity, diagonal values
% remove the LOI, for which ind_recurrence_points(i,1)==ind_recurrence_points(j,2), i=j
ptdiff = diff(ind_recurrence_points,1,2);
ind_recurrence_points_no_LOI = ind_recurrence_points(find(ptdiff),:);

if isempty(ind_recurrence_points_no_LOI)
    rqa_stat = zeros(1,7);
else
    %% Recurrence rate
    N = length(recurrence_plot);
    % REC = 1/(N^2)*sum(sum(recurrence_plot)); % with LOI
    REC = 1/(N^2-N)*size(ind_recurrence_points_no_LOI,1); % without LOI (Marwann2015)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal Line Structure Search
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diagonal lines represent such segments of the phase space trajectory
    % which run parallel for some time
    % Measures based on diagonal lines: DET, entropy, L, LMAX
    ind1 = 3;
    ind2 = 1;
    diagonal_lines_info = getLinePattern(ind_recurrence_points, ind1, ind2);
    
    % get the number of diagonals for diagonal size equal to or higher than
    % lmin:
    size_diagonal_lines = table2array(diagonal_lines_info(:,1));
    n_diagonal_lines = table2array(diagonal_lines_info(:,2));
    ind_lmin = find(size_diagonal_lines>=lmin);
    n_diagonal_lines_min = n_diagonal_lines(ind_lmin);
    size_diagonal_lines_min = size_diagonal_lines(ind_lmin);
    
    %% Entropy
    % reflects the complexity of the RP in respect of the diagonal lines,
    % e.g. for uncorrelated noise the value of ENTR is rather small,
    % indicating its low complexity
    
    prob = n_diagonal_lines/sum(n_diagonal_lines);
    nonz_prob = prob(find(prob));
    ENT = -sum(nonz_prob(ind_lmin).*(log2(nonz_prob(ind_lmin))));
    
    %% Determinism
    % ratio of recurrence points that form diagonal structures (of at least
    % length lmin) to all recurrence points
    % measure of predictability of the system
    DET = sum(size_diagonal_lines_min.*n_diagonal_lines_min)/sum(size_diagonal_lines.*n_diagonal_lines);
    
    %% Longest diagonal line found in the RP
    ind_LOI = size_diagonal_lines==N; % exclude LOI
    Lmax = max(size_diagonal_lines(~ind_LOI));
    
    %% Averaged diagonal line length
    ind_LOI_min = size_diagonal_lines_min==N; % exclude LOI
    L = sum(size_diagonal_lines_min(~ind_LOI_min).*n_diagonal_lines_min(~ind_LOI_min))/sum(n_diagonal_lines_min(~ind_LOI_min));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vertical Line Structure Search
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vertical lines represent segments which remain in the same phase
    % space region for some time.
    % Measures based on vertical lines: LAM, TT
    
    ind1 = 1;
    ind2 = 2;
    vertical_lines_info = getLinePattern(ind_recurrence_points, ind1, ind2);
    
    size_vertical_lines = table2array(vertical_lines_info(:,1));
    n_vertical_lines = table2array(vertical_lines_info(:,2));
    ind_lmin = find(size_vertical_lines>=lmin);
    n_vertical_lines_min = n_vertical_lines(ind_lmin);
    size_vertical_lines_min = size_vertical_lines(ind_lmin);
    
    
    %% Laminarity
    % ratio between the recurrence points forming the vertical structures
    % and the entire set of recurrence points (Marwan2007)
    LAM = sum(size_vertical_lines_min.*n_vertical_lines_min)/sum(size_vertical_lines.*n_vertical_lines);
    
    %% Trapping time
    % average length of vertical structures
    % how long the system remains in a specific state
    TT = sum(size_vertical_lines_min.*n_vertical_lines_min)/sum(n_vertical_lines_min);
    
    %%
    if isempty(REC)||isnan(REC); REC = 0; end
    if isempty(DET)||isnan(DET); DET = 0; end
    if isempty(Lmax)||isnan(Lmax); Lmax = 0; end
    if isempty(L)||isnan(L); L = 0; end
    if isempty(ENT)||isnan(ENT); ENT = 0; end
    if isempty(LAM)||isnan(LAM); LAM = 0; end
    if isempty(TT)||isnan(TT); TT = 0; end
    
    var_names = {'REC', 'DET', 'Lmax', 'L', 'ENT', 'LAM', 'TT'};

    rqa_stat = array2table(zeros(1,numel(var_names)), 'VariableNames', var_names);

    rqa_stat.Properties.Description = 'Recurrence Plot HRV Metrics';
    
    for ii = 1:numel(var_names)
        col_total_power = var_names{ii};
        eval(['rqa_stat{1,' num2str(ii) '} = ' var_names{ii} ';'])
        rqa_stat.Properties.VariableUnits{col_total_power} = '-';
        rqa_stat.Properties.VariableDescriptions{col_total_power} = var_names{ii};
    end
    
end
end

function line_size_frequency = getLinePattern(ind_recurrence_points, ...
    ind2inspect, ind2store)


%% (1) get the difference between the indices of the recurrence points
% draw a matrix of e.g. 10 by 10 and check that a diagonal line can be
% found by subtracting de indices of the rows to the columns, as the same
% value will be obtained for a diagonal line:
ind_recurrence_points = sortrows(ind_recurrence_points,1);
ind_recurrence_points = [ind_recurrence_points, ...
    ind_recurrence_points(:,2)-ind_recurrence_points(:,1)];

% sort the vector of differences to inspect each diagonal line which is
% given by sequences of values
ind_recurrence_points = sortrows(ind_recurrence_points,ind2inspect);


%% (2) Identify the number of points in each diagonal line. The row index
% of the points comprised in each diagonal line are stored in cell s


k = 1;
jj = 1;
for ii = 1:size(ind_recurrence_points,1)-1
    % store the row index of each point:
    
    % condition 1: the difference between row and column index of
    % recurrence point must be equal for two points
    % condition 2: the row indexes of the two points must be incremetal by
    % one
    s{k}(jj) = ind_recurrence_points(ii,ind2store);
    if ind_recurrence_points(ii,ind2inspect)==ind_recurrence_points(ii+1,ind2inspect) && ...
            ind_recurrence_points(ii+1,ind2store)==(ind_recurrence_points(ii,ind2store)+1)
        jj = jj+1;
    else
        jj = 1;
        k = k+1;
    end
    if ii==length(ind_recurrence_points)-1
        s{k}(jj) = ind_recurrence_points(ii+1,ind2store);
    end
    
end
s = s';


%% (3) count the size of each continuous line contained in each cell of s
% or, in other words, make the histogram of recurr_plot_no_LOI(:,ind2inspect)

size_lines = cellfun(@numel,s);

unique_line_size = unique(size_lines);

frequency_line_size = zeros(numel(unique_line_size),1);

for ii = 1:numel(unique_line_size)
    frequency_line_size(ii) = sum(size_lines==unique_line_size(ii));
end

line_size_frequency = table(unique_line_size, frequency_line_size, ...
    'VariableNames', {'Line size', 'Number of lines'});


end
