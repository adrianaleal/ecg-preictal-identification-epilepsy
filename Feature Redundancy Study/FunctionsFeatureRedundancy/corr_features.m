function [mat_corr,varargout] = corr_features(dataset,percentage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs: dataset sem missing values (dataset) e um threshold de separaçao (percentagem)
% Outputs: matrix de correlação (correlacao_features) e indices que estão acima do threshold (indices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the correlation of the features. 

% Inputs:
% dataset: feature dataset
% percentage: threshold of correlation

% Outputs:
% correlation_matrix
% indexes: indexes above the percentage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correlation_matrix = corrcoef(dataset, 'rows', 'pairwise');

correlacao = abs(correlation_matrix);

A = ones(size(dataset,2));
B = tril(A,-1);
mat_corr = correlacao(logical(B));


if nargin==2
    tt = tril(correlacao,-1); % removes the duplicate result of the correlation results between the windows 
    [x1,y1] = find(tt>=percentage);
    varargout = [x1,y1];
end

end
