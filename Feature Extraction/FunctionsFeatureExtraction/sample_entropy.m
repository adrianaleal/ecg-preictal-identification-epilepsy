function v_SampEn = sample_entropy(pv_DataIn,m,r_tolerance)


% ApEn reflects the likelihood that similar patterns of observations will
% not be followed by additional similar observations. (Wiki & http://physionet.incor.usp.br/physiotools/ApEn/)
% A time series containing many repetitive patterns has a relatively small
% ApEn; a less predictable process has a higher ApEn. (Wiki & http://physionet.incor.usp.br/physiotools/ApEn/)

% These are N raw data values from pv_DataIn equally spaced in time. (Wiki)
% HOWEVER, THERE IS NO TIME VECTOR INVOLVED IN THE COMPUTATION OF THE
% FEATURE... (adriana) (?????)


% Do we use instantaneous HR or the RR intervals ????
% Physionet uses the HR (http://physionet.incor.usp.br/physiotools/ApEn/)

% r defines the criterion of similarity.

% WEAKNESSES (http://physionet.incor.usp.br/physiotools/ApEn/):
% (i) strong dependence on sequence length
% (ii) poor self-consistency (i.e., the observation that ApEn for one data
% set is larger than ApEn for another for a given choice of m and r should,
% but does not, hold true for other choices of m and r).

if (m < 0 || r_tolerance < 0)
    error('Invalid parameter values');
end

[N, ~] = size(pv_DataIn);

if N == 1
    pv_DataIn = pv_DataIn';
    N = length(pv_DataIn);
end

rs = std(pv_DataIn);
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
    % P = 0; % for ApEn
    
    % Initialize template-match counters. A is the number of template matches
    % of length m+1, and B is the number of template matches of length m.
    A = 0;
    B = 0;
%     A = zeros(N-m,1);
%     B = A;
    % count_j = [];
    % Bm(r) is then the probability that two sequences will match for m 
    % points, whereas Am(r) is the probability that two sequences will 
    % match for m+1 points. (Richman2000)

    
    % SampEn does not use a templatewise approach, and A and B cumulate
    % for all the templates.
    for ii = 1:N-m
%         Bi = 0;
%         Ai = 0;
        for jj = (ii+1):N-m+1
            % count_j = [count_j; ii jj];
            D1 = abs(pv_DataIn(ii)-pv_DataIn(jj));
            D2 = abs(pv_DataIn(ii+1)-pv_DataIn(jj+1));
            
            if(D1>D2)
                DA = D1;
            else
                DA = D2;
            end
            
            if DA<=r
                B = B+1;
                if jj<=N-m
                    D3 = abs(pv_DataIn(ii+m)-pv_DataIn(jj+m));
                    if D3<=r
                        A = A+1;
                    end
                end
            end
        end
        
        % B(ii) = Bi;
        %bA(ii) = Ai;
        % for ApEn:
        % if (A>0)&&(B>0)
        %     Pi = (B/A);
        %     P = P + log(Pi);
        % end
    end
    
%     A
%     B
    %     ApEn =(-2)*P*(1/(N-m));
    
    v_SampEn = - log((A/B)*((N-m+2)/(N-m)));
    
    % A = A(r)*(N-(m+1)-1)*(N-(m+1))/2
    % B = B(r)*(N-m-1)*(N-m)/2
    
    %     disp(['Original Approximate Entropy = ' num2str(ApEntropy)])
    %
    %     disp('****************************************************************')
end

end
















% WHEN COMPARING TO FUNCTION sampen.m
%
% matches = NaN(m+1,N);
% for ii = 1:1:m+1
%     matches(ii,1:N+1-ii) = pv_DataIn(ii:end);
% end
% matches = matches';
% d_m = pdist(matches(:,1:m), 'Chebychev')';
% d_m1 = pdist(matches(:,1:m+1), 'Chebychev')';
%
% matches2 = NaN(N,m+1);
% d_m2 = [];
% d_m22 = [];
% for ii = 1:N
%     matches2(ii,1) = pv_DataIn(ii);
%     if ii+1<=N
%         matches2(ii,2) = pv_DataIn(ii+1);
%         if ii+2<=N
%             matches2(ii,3) = pv_DataIn(ii+2);
%             for jj = (ii+1):N-m+1
%                 D1 = abs(pv_DataIn(ii)-pv_DataIn(jj));
%                 D2 = abs(pv_DataIn(ii+1)-pv_DataIn(jj+1));
%
%                 if(D1>D2)
%                     DA = D1;
%                 else
%                     DA = D2;
%                 end
%
%                 d_m2 = [d_m2; DA];
%
%                 if jj+2<=N
%                     D3 = abs(pv_DataIn(ii+2)-pv_DataIn(jj+2));
%                     if DA<D3
%                         DA = D3;
%                     end
%                     d_m22 = [d_m22; DA];
%                 end
%
%             end
%
%
%         end
%     end
% end
