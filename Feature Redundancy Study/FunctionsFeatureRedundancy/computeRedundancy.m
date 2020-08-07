function computeRedundancy(n_feat, n_comb, dim, variables2compute, ...
    dimension_string, feature_dataset, comb, chosenDim, ...
    n_variables2compute, subfolder2save)


normalityAssessment = zeros(n_feat, dim);
pNormalityAssessment = normalityAssessment;


% initialize variables2compute
eval([variables2compute{1} ' = zeros(n_comb, dim);'])
for ll = 1:n_variables2compute-1
    eval([variables2compute{ll} ' = ' variables2compute{1} ';'])
end

for ii = 1:dim
    
    disp([dimension_string ' ' num2str(ii)])
    
    %% Compute correlation between features (Pearson's)
    
    if strcmp(chosenDim, 'seizures')
        data2compare = squeeze(feature_dataset(ii,:,:))';
    elseif strcmp(chosenDim, 'windows')
        data2compare = squeeze(feature_dataset(:,:,ii));
    end
    
    % get the Pearson's correlatio coefficient for each dimension:
    if any(strcmp(variables2compute, 'correlation_between_features_pearson'))
        correlation_between_features_pearson(:,ii) = corr_features(data2compare);
    end
    
    %% Check feature normality
    for ff = 1:n_feat
        [normalityAssessment(ff,ii), pNormalityAssessment(ff,ii)] = kstest(data2compare(:,ff));
    end
    % returns a test decision for the null hypothesis that the data in vector
    % x comes from a standard normal distribution, against the alternative that
    % it does not come from such a distribution, using the one-sample
    % Kolmogorov-Smirnov test. The result h is 1 if the test rejects the null
    % hypothesis at the 5% significance level, or 0 otherwise.
    
    
    %% Compute correlation and mutual information between features
    
    for jj = 1:n_comb
        getMat2correlate = data2compare(:,comb(jj,:));
        getMat2correlateNoNaN = getMat2correlate(~sum(isnan(getMat2correlate),2)>0,:);
        
        % Spearman's correlation **************************************
        if any(strcmp(variables2compute, 'correlation_between_features_spearman'))
            [RHO,PVAL] = corr(getMat2correlateNoNaN(:,1), getMat2correlateNoNaN(:,2), ...
                'Type','Spearman');
            correlation_between_features_spearman(jj,ii) = RHO;
        end
        
        % Kendall's correlation ***************************************
        if any(strcmp(variables2compute, 'correlation_between_features_kendall'))
            [RHO,PVAL] = corr(getMat2correlateNoNaN(:,1), getMat2correlateNoNaN(:,2), ...
                'Type','Kendall');
            correlation_between_features_kendall(jj,ii) = RHO;
        end
        
        % Average mutual information **********************************
        if any(strcmp(variables2compute, 'mutual_information_between_features_AMI'))
            plotFigure = 0;
            mutual_information_between_features_AMI(jj,ii) = average_mutual_information( ...
                getMat2correlateNoNaN, plotFigure);
        end
        
    end
    
end

for ll = 1:n_variables2compute
    eval(['save(fullfile(subfolder2save, ''' variables2compute{ll} ...
        '.mat''), ''' variables2compute{ll} ''')'])
end

end