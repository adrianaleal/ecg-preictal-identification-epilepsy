function [features_segment, comp_time_segment] = HRV_features(RR_interval_signal_segment, ...
    time_RR_intervals, features_segment, comp_time_segment, ll, ...
    segment_size_seconds, plotFigure)

if plotFigure
    pat_seiz_name = 'pat_8902_seiz_1_win_1940';
end

% countAbnormalRR = features_segment.nAbnormalRR(ll);
% if (countAbnormalRR/numel(RR_interval_signal_segment))>0.8
%     return
% end

% *****************************************************************
% to compute the classical scalar HRV parameters (SDNN, pNN50, RMSSD)
% there is no need for the instantaneous heart rate, only inter-beat RR
% intervals

    
%% Time-domain features

% disp('>> Computing time-domain features')


if segment_size_seconds>=120
    tic
    features_segment = time_features(RR_interval_signal_segment,...
        features_segment, ll);
    
    % fprintf('[%.3f sec] >> Elapsed time for time-domain features \n', toc);
%     fprintf('\n')
    comp_time_segment.time(ll) = toc;
    
end



%% Frequency-domain features

% disp('>> Computing frequency-domain features')

% [hrv_fd, pxx, f_axis, plot_data] = hrv_freq(RR_interval_signal_segment);


if segment_size_seconds>=120%size(RR_interval_signal_segment,2)>17
    tic
    % plot_frequency = 1;
    
    [hrv_fd, interp_RRI_signal_segment] = frequency_features(RR_interval_signal_segment,...
        time_RR_intervals,segment_size_seconds, plotFigure);
    
    if ~isempty(hrv_fd)
        % Get the relative power in each band for the lomb method
        
        features_segment.TOTAL_POWER(ll) = hrv_fd{1,1};
        features_segment.VLF_POWER(ll) = hrv_fd{1,2};
        features_segment.LF_POWER(ll) = hrv_fd{1,3};
        features_segment.HF_POWER(ll) = hrv_fd{1,4};
        features_segment.VLF_NORM(ll) = hrv_fd{1,5};
        features_segment.LF_NORM(ll) = hrv_fd{1,6};
        features_segment.HF_NORM(ll) = hrv_fd{1,7};
        features_segment.LFtoHF(ll) = hrv_fd{1,8};
    end

    comp_time_segment.frequency(ll) = toc;
    
    if plotFigure
        saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
            [pat_seiz_name '_freq']))
        close all
    end
    
else
    interp_RRI_signal_segment = [];
end

% fprintf('[%.3f sec] >> Elapsed time for frequency-domain features \n', toc);
% fprintf('\n')


%% Nonlinear features

% disp('>> Computing nonlinear features')

% fprintf('\n')


% Poincaré plot

% disp('>> Computing Poincaré plot')
% - SD1: Poincare plot SD1 descriptor (std. dev. of intervals along
% the line perpendicular to the line of identity).
% - SD2: Poincare plot SD2 descriptor (std. dev. of intervals along
% the line of identity).
% tic
% plot_poincare = 1;


if segment_size_seconds>=120
    tic
    [SD1, SD2, SDRR] = poincare_plot(RR_interval_signal_segment,plotFigure);
    features_segment.SD1(ll) = SD1;
    features_segment.SD2(ll) = SD2;
    features_segment.SD1toSD2(ll) = SD1/SD2;
    comp_time_segment.poincare(ll) = toc;
    
    if plotFigure
        saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
            [pat_seiz_name '_poincare']))
        close all
    end
    
end

% fprintf('[%.3f sec] >> Elapsed time for Poincaré plot \n', toc);
% fprintf('\n')



%%
% disp('>> Compute Detrended Fluctuation Analysis')
% - alpha1: Log-log slope of DFA in the low-scale region.
% - alpha2: Log-log slope of DFA in the high-scale region.
% Limitations: at least 8000 data points should be used; monofractal method;
% normal-to-normal interbeat intervals are required; dependency on editing
% ectopic beats.Voss2009

n_beats_minimum = 100;
if size(RR_interval_signal_segment,2)>=(n_beats_minimum-1)
    tic
    %         tnn = [RR_interval_signal_segment(1) cumsum(RR_interval_signal_segment(1:end-1))]/1e3'; %
    [~, ~, alpha1, alpha2] = detrended_fluctuation_analysis(time_RR_intervals, ...
        RR_interval_signal_segment, plotFigure);
    features_segment.DFA_alpha1(ll) = alpha1;
    features_segment.DFA_alpha2(ll) = alpha2;
    comp_time_segment.DFA(ll) = toc;
    
    if plotFigure
        saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
            [pat_seiz_name '_dfa']))
        close all
    end
end
% fprintf('[%.3f sec] >> Elapsed time for detrended fluctuation analysis \n', toc);
% fprintf('\n')

%%
% tic
% disp('>> Computing approximate entropy')
%
% n_beats_min = 200;% Shaffer2017, Kuusela2012
% if size(RR_interval_signal_segment,2)>=(n_beats_min-1)
%     ApEn = approximate_entropy(RR_interval_signal_segment,2,0.2);
%     features_segment.ApEn(ll) = ApEn;
% end
% fprintf('[%.3f sec] >> Elapsed time for approximate entropy \n', toc);


%%
% disp('>> Compute multiscale entropy')
% Which means computing the sample entropy for the mean of
% differently sized windows (different scales)
% May be feasible for an entire window which can be divided in a
% minimum of 100 beats windows...

% n_beats_minimum = 100;
% if size(RR_interval_signal_segment,2)>=(n_beats_minimum-1)
%     tic
%     m = 2;
%     r_tolerance = 0.2;
% %     mse_max_scale = round(numel(RR_interval_signal_segment)/n_beats_minimum);
% %     [mse_result, ~] = multi_scale_entropy(RR_interval_signal_segment,...
% %         'mse_max_scale', mse_max_scale, 'sampen_r', r_tolerance, ...
% %         'sampen_m', m, 'plot', plotFigure);
% %     SampEn = mse_result(1);
%     
% %     tic
%     SampEn = sample_entropy(RR_interval_signal_segment, m, r_tolerance);
% %     toc
% %     tic
% %     SampEn = sampen(RR_interval_signal_segment, m, r_tolerance);
% %     toc
%     
%     features_segment.SampEn(ll) = SampEn;
%     comp_time_segment.SampEn(ll) = toc;
% end
% % fprintf('[%.3f sec] >> Elapsed time for multiscale entropy \n', toc);
% % fprintf('\n')



%% normalize signal

if ~isempty(interp_RRI_signal_segment)
    
    
    sig2analyse = interp_RRI_signal_segment; % RR_interval_signal_segment
    % sig_norm = (sig2analyse + ...
    %     abs(min(sig2analyse)))/max(sig2analyse + ...
    %     abs(min(sig2analyse))); %(read Epileptic Seizure and the EEG, Varsavsky2011, page 137)
    
    sig_norm = detrend(sig2analyse)/std(sig2analyse); % with this an error
    % is obtained in correlation dimension
    
    
    %%
    
    n_beats_min = 200;% Shaffer2017
    if size(RR_interval_signal_segment,2)>=(n_beats_min-1)
        tic
        % estimates the approximate entropy of the uniformly sampled (matlab) 
        % time-domain signal X by reconstructing the phase space
        % The default value of Radius is, 0.2*variance(X), if X has a 
        % single column.
        % m default value is 2
        approxEnt = approximateEntropy(sig_norm);
        
        features_segment.ApEn(ll) = approxEnt;
        comp_time_segment.ApEn(ll) = toc;
        
        % the results are similar to the function bellow:
%         ApEn = approximate_entropy(sig_norm,2,0.2);
        
    end
    
    n_beats_minimum = 100;
    if size(RR_interval_signal_segment,2)>=(n_beats_minimum-1)
        %% disp('>> Compute Sample Entropy')
        tic
        m = 2;
        r_tolerance = 0.2;
        SampEn = sample_entropy(sig_norm, m, r_tolerance);
        features_segment.SampEn(ll) = SampEn;
        comp_time_segment.SampEn(ll) = toc;
        
        
        %%
%         disp('>> Compute Correlation Dimension and Largest Lyapunov Exponent ')
        % tic
        
        %% compute time delay based on the autocorrelation
        % tau = getEmbeddingTimeDelay(sig_norm, 'mutual information');
        
        %% Embedding dimension, calculated with False Nearest Neighbors(FNN)
        
        % RPS = RecPhaseSpace(sig_norm,tau);
        % [eDim, ed, E] = CalcMinEmbDimension(sig_norm, RPS, tau, 'euclidean', 1);
        
        
        % phaseSpaceReconstruction requires that sig_norm is uniformly sampled
        tic
        [attractorMat, tauMat, eDimMat] = phaseSpaceReconstruction(sig_norm);
        features_segment.tau_rec_phase_space(ll) = tauMat;
        features_segment.embDim_rec_phase_space(ll) = eDimMat;
        PSR_comp_time = toc;
        
        if plotFigure
            figure()
            phaseSpaceReconstruction(sig_norm, tauMat, eDimMat);
            close all
        end
        
        % RPS = RecPhaseSpace(sig_norm,tauMat);
        % ed = CalcDistBPoints(RPS, tauMat);
        % E = CalcFalseNeighbours(RPS,tauMat,ed);
        % [largLyapunovExp, attractor] = CalcLyapExponent(sig_norm, ed, eDimMat, ...
        %     tauMat, RPS, 0);
        
        % The Lyapunov exponent measures the rate at which the subsequent orbits
        % diverge, thus measuring the amount of stability or instability in the
        % system or its sensitivity to initial conditions.
        % [largLyapunovExp2, corrDim] = lyap_compu(sig_norm, tauMat, eDimMat);
        
        
        
        % [ ns, lnd ] = lyapunov(sig_norm, eDimMat, tauMat, 357, 5, 10, 10);
        % d = exp(1).^lnd
        
        % INPUT:
        %   x -- scalar time series
        %   d -- embedding dimension
        %   m -- delay parameter
        %   N -- number of initial (fiducial) points to use (should be large)
        %   R -- number of nearest neighbors to use (R = 1 should work)
        %   n -- maximum number of steps ahead to proceeed
        %   p -- (optional) supply 0 to not remove temporal correlations
        %        or number of points to remove. The default is p=d*m.
        %        NOTE: only need to remove these correlations for flow data!!
        tic
        lyapExpMat = lyapunovExponent(sig_norm, 4, 'Dimension', eDimMat, ...
            'Lag', tauMat, 'ExpansionRange', 10, 'MinSeparation', 125);
        features_segment.largLyapunov(ll) = lyapExpMat;
        comp_time_segment.largLyapunov(ll) = toc+PSR_comp_time;
        
        
        if plotFigure
            if eDimMat==2
                figure()
                plot(attractorMat(:,1), attractorMat(:,2))
                xlabel('$x(k)$')
                ylabel(['$x(k-' num2str(tauMat) ')$'])
            elseif eDimMat==3
                % A three-dimensional embedding uses a delay of 2*tau for the third
                % dimension.
                figure()
                plot3(attractorMat(:,1), attractorMat(:,2), attractorMat(:,3))
                xlabel('$x(k)$')
                ylabel('$x(k-\tau)$')
                zlabel('$x(k-2\tau)$')
            end
            
            figure()
            lyapunovExponent(sig_norm, 4, 'Dimension', eDimMat, ...
                'Lag', tauMat, 'ExpansionRange', 10, 'MinSeparation', 125);

            saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
                [pat_seiz_name '_largLyapunovExp']))
            close all
        end
        
        
        %% estimates the correlation dimension of the uniformly sampled 
        % time-domain signal
        
        if plotFigure
            correlationDimension(sig_norm,'Dimension', eDimMat, 'Lag',tauMat, ...
                'MaxRadius', 0.05, 'NumPoints', length(sig_norm)); % ,'NumPoints',Np
        end
        tic
        [corrDimMat,rRange,corInt] = correlationDimension(sig_norm, ...
            'Dimension', eDimMat, 'Lag', tauMat, ...
            'NumPoints', length(sig_norm));
        features_segment.corrDim(ll) = corrDimMat;
        comp_time_segment.corrDim(ll) = toc+PSR_comp_time;
        
        if plotFigure
            saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
                [pat_seiz_name '_corrDim']))
            close all
        end
    
        % fprintf('[%.3f sec] >> Elapsed time for Correlation Dimension and Largest Lyapunov Exponent \n', toc);
        fprintf('\n')
        
        %%
%         disp('>> Compute Recurrence Quantification Analysis')
        tic
        [rqa_stat, ~, recurdata] = recurrenceAnalysis(sig_norm, tauMat, ...
            eDimMat, attractorMat, plotFigure);
        
        % rqa_stat = [REC DET Lmax L ENT LAM TT];
        
        features_segment.RQA_REC(ll) = rqa_stat(1);
        features_segment.RQA_DET(ll) = rqa_stat(2);
        features_segment.RQA_Lmax(ll) = rqa_stat(3);
        features_segment.RQA_L(ll) = rqa_stat(4);
        features_segment.RQA_ENT(ll) = rqa_stat(5);
        features_segment.RQA_LAM(ll) = rqa_stat(6);
        features_segment.RQA_TT(ll) = rqa_stat(7);
        comp_time_segment.RQA(ll) = toc+PSR_comp_time;
        % fprintf('[%.3f sec] >> Elapsed time for recurrence quantification analysis \n', toc);
        fprintf('\n')
        
        if plotFigure
            eval([pat_seiz_name '_rqa_distance_mat = recurdata;'])
            save(fullfile(cd, 'featureExtractionImages', ...
                [pat_seiz_name '_rqa_distance_mat.mat']))
            
            saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
                [pat_seiz_name '_rqa_RP']))
            
        end
    end
    
else
    disp('Small RR signal!!')
end

fprintf('\n')
fprintf('\n')


end
