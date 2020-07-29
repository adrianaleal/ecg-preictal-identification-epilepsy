function firstLocalMinAMI = getAverageMutualInformation(sig, maxTau, ...
    plotFigure, method)

if nargin<4
    %% COMPARE ALL METHODS TO OBTAIN AVERAGE MUTUAL INFORMATION
    
    %% amutual from tstoolbox/@signal/amutual
    s = signal(sig);
    AMI_obj = amutual(s, maxTau);
    % amutual(signal, maxtau, bins) | maxtau = 32 and number of bins 128
    AMI = data(AMI_obj);
    [firstLocalMinAMI, AMIvec] = findLocalMinimum(AMI); % , AMI_obj
    
    %% 
    nbins = 128;% default in amutual
    AMI2 = ami(sig,0,maxTau,nbins);
    [firstLocalMinAMI2, AMI2vec] = findLocalMinimum(AMI2);

    %%
    % AMI3: "An efficient algorithm for the computation of average mutual 
    % information: Validation and implementation in Matlab" by Thomas 2014
    
    % get the matrix of delayed versions of sig
    nmax = size(sig,1) - maxTau; % <- that's the maximum size possible!
    n = nmax;
    AMI3 = zeros(maxTau+1,1);
    AMI4 = zeros(maxTau+1,1);
    for ii = 0:maxTau
        if ii>0
            AMI3(ii) = average_mutual_information([no_delay_vec, sig(1+ii:n+ii)]);
            [v,~] = ami_v2(no_delay_vec,sig(1+ii:n+ii),ii);
            AMI4(ii) = v;
        else
            no_delay_vec = sig(1+ii:n+ii);
        end
    end

    [firstLocalMinAMI3, AMI3vec] = findLocalMinimum(AMI3);
    [firstLocalMinAMI4, AMI4vec] = findLocalMinimum(AMI4);
    
    %%
    
    figure()
    plot(AMI), hold on
    plot(AMI2)
    plot(AMI3)
    plot(AMI4)
    plot(AMIvec, AMI(AMIvec),'k*')
    plot(AMI2vec, AMI2(AMI2vec),'k*')
    plot(AMI3vec, AMI3(AMI3vec),'k*')
    plot(AMI4vec, AMI4(AMI4vec),'k*')
    ylabel('Mutual information')
    xlabel('$\tau$')
    axis tight
    hold off
    legend('AMI tstoolbox', 'AMI2', 'AMI Thomas2014', 'AMI Grinsted2006', 'Local minimums')


else
    
    if strcmp(method, 'tstoolbox')
        % amutual from tstoolbox/@signal/amutual: Auto mutual information ????
        % Mutual information between a signal and its time-delayed version.
        
        s = signal(sig);
        AMI_obj = amutual(s, maxTau);
        % amutual(signal, maxtau, bins) | maxtau = 32 and number of bins 128
        AMI = data(AMI_obj);
        
    elseif strcmp(method, 'AMI2')
        nbins = 128;% default in amutual
        AMI = ami(sig,0,maxTau,nbins);
        
    else
        
        if strcmp(method, 'Thomas2014')
            AMI = zeros(maxTau+1,1);
        elseif strcmp(method, 'Grinsted2006')
            AMI = zeros(maxTau+1,1);
        end
        
        nmax = size(sig,1) - maxTau; % <- that's the maximum size possible!
        n = nmax;
        for ii = 1:maxTau
            no_delay_vec = sig(1:end-ii); % instead of no_delay_vec = sig(1:n);
            delay_vec = sig(1+ii:end); % instead of delay_vec = sig(1+ii:n+ii);
            if strcmp(method, 'Thomas2014')
                AMI(ii) = average_mutual_information([no_delay_vec, delay_vec]);
            elseif strcmp(method, 'Grinsted2006')
                [v,~] = ami_v2(no_delay_vec, delay_vec, ii);
                AMI(ii) = v;
            end
        end
    end
    
    % to compare different methods of search for local minimuns add 
    % ', AMI_obj' to the input arguments:
    [firstLocalMinAMI, AMIvec] = findLocalMinimum(AMI);
    
    
    if plotFigure   
        figure()
        plot(AMI), hold on
        plot(AMIvec, AMI(AMIvec),'k*'), hold off
        ylabel('Average mutual information')
        xlabel('$\tau$')
        axis tight
        legend(['AMI ' method], 'Local minimums')
    end
    
end


end
