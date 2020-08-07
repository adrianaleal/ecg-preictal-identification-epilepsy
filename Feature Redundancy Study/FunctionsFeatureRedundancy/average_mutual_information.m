function AMI = average_mutual_information(data, plotFigure)
% function AMI = average_mutual_information(data)
% Calculates average mutual information between two columns of data. It
% uses kernel density estimation, with a globally adjusted Gaussian kernel.
%
% Input: n-by-2 matrix, with data sets in adjacent column vectors.
%
% Output: scalar.

% Source: "An efficient algorithm for the computation of average mutual
% information: Validation and implementation in Matlab" by Thomas 2014


n = length(data);
X = data(:,1);
Y = data(:,2);

% Example below is for normal reference rule in 2 dims, Scott (1992).
hx = std(X)/(n^(1/6));
hy = std(Y)/(n^(1/6));

% Estimated univariate marginal probability density functions:
P_x = univariate_kernel_density(X, X, hx, plotFigure);
P_y = univariate_kernel_density(Y, Y, hy, plotFigure);

% Estimated joint probability density function:
JointP_xy = bivariate_kernel_density(data, data, hx, hy, plotFigure);

% If the log base 2 is used, the units of mutual information are bits.
AMI = sum(log2(JointP_xy./(P_x.*P_y)))/n;

% IMPORTANT: it is not possible to obtain entropies from the marginal
% probability density functions see second paragraph, first column, page
% two in Thomas2014

end

function y = univariate_kernel_density(value, data, window, plotFigure)

% function y = univariate_kernel_density(value, data, window)
% Estimates univariate density using kernel density estimation.

% Inputs:
% value (m-vector), where density is estimated;
% data (n-vector), the data used to estimate the density;
% window (scalar), used for the width of density estimation.

% Output: m-vector of probabilities.

h = window;
n = length(data);
m = length(value);

% We use matrix operations to speed up computation of a double-sum.

% Prob = zeros(n, m);
G = Extended(value, n); % column copies of value
H = Extended(data', m); % row copies of data

input2Dmat = (G - H)/h;

 
Prob = normpdf(input2Dmat);
% Default values for MU and SIGMA are 0 and 1 respectively.

% see univariate distributions (in each column)
if plotFigure
    figure()
    plot(input2Dmat,Prob)
end


fhat = sum(Prob)/(n*h);

y = fhat';

end

function y = bivariate_kernel_density(value, data, Hone, Htwo, plotFigure)
% function y = bivariate_kernel_density(value, data, Hone, Htwo)
% Calculates bivariate kernel density estimates of probability.

% Inputs:
% value (m x 2 matrix), where density is estimated;
% data (n x 2 matrix), the data used to estimate the density;
% Hone (scalar) and Htwo (scalar) to use for the widths of density estimation.

% Output: m-vector of probabilities estimated at the values in 'value'.

s = size(data);
n = s(1);
t = size(value);
number_pts = t(1);

% rho_matrix = corr(data); % original code
rho_matrix = corrcoef(data);
% another option with the same result:
% rho_matrix = corrcoef(data,'rows','complete');

rho = rho_matrix(1,2);

% The adjusted covariance matrix:
W = [Hone^2 rho*Hone*Htwo; rho*Hone*Htwo Htwo^2];

Differences = linear_depth(value,-data);
% [D, PD] = allfitdist(Differences(:,1), 'PDF');
% [D, PD] = allfitdist(Differences(:,2), 'PDF');

% Multivariate normal probability density function:
try
    prob = mvnpdf(Differences,[0 0],W); % original code
catch
    % ERROR: SIGMA (W) must be a square, symmetric, positive definite matrix.
    % Meaning: this happens if the diagonal values of the covariance matrix
    % are (very close to) zero.
    W2 = W + .0001 * eye(2);
    prob = mvnpdf(Differences,[0 0],W2);
    
    % check if W is a symmetric positive definite matrix
    % d = eig(W)
    % isposdef = all(d) > 0
end

Cumprob = cumsum(prob);


if plotFigure
    
    % Plot the probability density values.
    figure
    scatter3(Differences(:,1),Differences(:,2),prob)
    xlabel('X1')
    ylabel('X2')
    zlabel('Probability Density')
    
end


y = zeros(number_pts,1);
y(1) = (1/n)*Cumprob(n);
for ii = 2:number_pts
    index = n*ii;
    y(ii) = (1/(n))* (Cumprob(index)-Cumprob(index - n));
end

end

function y = linear_depth(feet, toes)
% linear_depth takes a matrix 'feet' and lengthens it in blocks, takes a
% matrix 'toes' and lengthens it in Extended repeats, and then adds the
% lengthened 'feet' and 'toes' matrices to achieve all sum combinations of
% their rows.
% feet and toes have the same number of columns

if size(feet, 2) == size(toes, 2)
    a = size(feet, 1);
    b = size(toes, 1);
    Blocks = zeros(a*b, size(toes, 2));
    Bricks = Blocks;
    
    for ii = 1:a
        Blocks((ii-1)*b + 1: ii*b,:) = Extended(feet(ii,:),b);
        Bricks((ii-1)*b + 1: ii*b,:) = toes;
    end
end

y = Blocks + Bricks;

end

function y = Extended(vector,n)
% Takes an m-dimensional row vector and outputs an n-by-m matrix with
% n-many consecutive repeats of the vector. Similarly, it takes an
% m-dimensional column vector and outputs an m-by-n matrix.
% Else, it returns the original input.

M = vector;

if size(vector,1) == 1
    M = zeros(n,length(vector));
    for ii = 1:n
        M(ii,:) = vector;
    end
end

if size(vector,2) == 1
    M = zeros(length(vector),n);
    for ii = 1:n
        M(:,ii) = vector;
    end
end

y = M;

end