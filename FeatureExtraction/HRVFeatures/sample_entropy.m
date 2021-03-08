function SampEn = sample_entropy(input_signal, m, r_tolerance)

% Compute sample entropy

% Inputs:
% - input_signal (double): RR interval series
% - m (double): 
% - r_tolerance (double): defines the criterion of similarity

% Outputs:
% - SampEn: sample entropy

% READ:
% (1) Richman2000: Richman, J. S. & Moorman, J. R. Physiological 
% time-series analysis using approximate entropy and sample entropy. Am. J.
% Physiol. Circ. Physiol. 278, H2039–H2049, DOI: 
% 10.1152/ajpheart.2000.278.6.H2039 (2000).
% (2) Aboy2007: Aboy, M., Cuesta-Frau, D., Austin, D. & Mico-Tormos, P. 
% Characterization of Sample Entropy in the Context of Biomedical Signal 
% Analysis. In 2007 29th Annual International Conference of the IEEE 
% Engineering in Medicine and Biology Society, 5942–5945, DOI: 
% 10.1109/IEMBS.2007.4353701 (IEEE, 2007).


%%
if (m < 0 || r_tolerance < 0)
    error('Invalid parameter values');
end

[N, ~] = size(input_signal);

if N == 1
    input_signal = input_signal';
    N = length(input_signal);
end

rs = std(input_signal);
if ~isempty(rs)
    r = r_tolerance*rs; % r specifies a filtering level or tolerance for accepting matches.
end
% It is convenient to set the tolerance as r*SD, the standard deviation
% of the data set, allowing measurements on data sets with different
% amplitudes to be compared (Richman2000).
% If r is too small, noise affects the SampEn measure. If r is too
% large, some changes of the signal are not detected (Aboy2007).

%% Two ways of computing the approximate entropy:

%% By not considering repeated matching and self-matching
% And lacking information about two samples because the same size of
% the distance matrices is needed to compare them
% Original computation

if m==2
    
    % Initialize template-match counters. A is the number of template matches
    % of length m+1, and B is the number of template matches of length m.
    A = 0;
    B = 0;
    % Bm(r) is then the probability that two sequences will match for m 
    % points, whereas Am(r) is the probability that two sequences will 
    % match for m+1 points. (Richman2000)

    
    % SampEn does not use a templatewise approach, and A and B cumulate
    % for all the templates.
    for ii = 1:N-m

        for jj = (ii+1):N-m+1
            D1 = abs(input_signal(ii)-input_signal(jj));
            D2 = abs(input_signal(ii+1)-input_signal(jj+1));
            
            if(D1>D2)
                DA = D1;
            else
                DA = D2;
            end
            
            if DA<=r
                B = B+1;
                if jj<=N-m
                    D3 = abs(input_signal(ii+m)-input_signal(jj+m));
                    if D3<=r
                        A = A+1;
                    end
                end
            end
        end
    end
 
    SampEn = - log((A/B)*((N-m+2)/(N-m)));
end

end
