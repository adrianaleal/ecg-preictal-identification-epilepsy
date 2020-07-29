function [lyap_r, corr_dim] = lyap_compu(sig, tau, embed_dim)


if size(sig,2)~=1
    sig = sig';
end

s = signal(sig);

if nargin<3
%% Cao's algorithm (to find embedding dimension)
max_dim = 10;

c = cao(s,max_dim,tau,5,-1); %cao(signal, maxdim, tau, NNR, Nref)
% | maxdim accepted = 10 | tau = first minimum of auto mutual information
% method | NNR = number of nearest neighbors = 5 | Nref - number of 
% reference points (-1 means: use all points)

% derivative of cao's method plot (if under 0.1 => x = embedding dimension)
embed_dim = diff(data(c),1);
embed_dim = reverse(signal(embed_dim));
offset = 0.1;
embed_dim = minus(embed_dim,offset); %returns the 'kink in the cao's graph'

% Define embedding dimension according to Cao's method
if ((min(embed_dim)<0) && (max(embed_dim)>0))
    f_zero = round(firstzero(embed_dim));
    embed_dim = (max_dim-f_zero);
else
    f_zero = max_dim;
    embed_dim = 6;
end

fprintf('Embedding Dimension: %d\n', embed_dim)
end

%% e = embbeded time series | attractor
e = embed(s, embed_dim, tau);
% embed(signal, dimension, delay|tau, windowtype)

%% Lyapunov exponent
l = largelyap(e,-1,10,10,5);
% largelyap (signal, number of randomly chosen reference points (-1 means: 
% use all points), stepsahead | maximal length of prediction in samples, 
% past | number of points excluded in same trajectory, nnr | number of 
% nearest neighbours)


% largest lyapunov exponent | slope of linear segment of the largelyap plot
[x_l, y_l] = size(data(l));

counter = 0;
lyap_data = data (l);

error_time = 10;

for jj = 2:x_l
    if (lyap_data(jj) - lyap_data(jj-1) < 0.05)
        counter = counter + 1;
    else
        counter = 0;
    end
    
    if (counter == 3)
        error_time = jj;
        break;
    end
end

%  l = slope of the linear section
lyap_r = lyap_data(error_time)/error_time;


%% Correlation dimension
corr_dim = takens_estimator(e,-1,0.05,40);
% n - number of randomly chosen reference points (n == -1 means: use all points)
% range - maximal relative search radius (relative to attractor size) 0..1
% past - number of samples to exclude before and after each reference index


end