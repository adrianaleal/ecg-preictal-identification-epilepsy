function [rqa_stat, attractor, recurdata] = recurrenceAnalysis(sig, tau, ...
    eDim, attractor, plotFigure)

% % mutual information test to determine the time delay
% mi = mutual(sig,[], [], 1);
% 
% % fnn test to determine the embedding dimension
% maxDim = 10;
% out = false_nearest(sig,1,maxDim,tau);
% 
% 
% fnn = out(:,1:2);
% figure();
% plt = plot(fnn(:,1),fnn(:,2),'o-','MarkerSize',4.5);
% title('False nearest neighbor test');
% xlabel('dimension');
% ylabel('FNN');
% grid on;

% phase space plot
if nargin<4
    attractor = phasespace(sig, tau, eDim, 1); % tau = 3, eDim = 8
end
% if size(attractor,2)<4
% figure();
% if size(attractor,2)==3
%     plot3(attractor(:,1),attractor(:,2),attractor(:,3),'-');
%     zlabel('x(t+2$\tau$)');
% elseif size(attractor,2)==2
%     plot(attractor(:,1),attractor(:,2),'-');
% end
% title('EKG time-delay embedding - state space plot');
% grid on;
% xlabel('x(t)');
% ylabel('x(t+$\tau$)');
% end


% color recurrence plot
recurdata = cerecurr_y(attractor, plotFigure);% plotFigure


% epsilon = sqrt(eDim)*
% black-white recurrence plot

% In: http://www.recurrence-plot.tk/forum/viewtopic.php?f=1&t=3583&p=5575#p5575

% maximal phase space diameter = |xMax – xMin|:
maxDiam = abs(max(attractor(:))-min(attractor(:)));
% maxDiam = pss(x,m,tau); in http://forum.recurrence-plot.tk/viewtopic.php?f=1&t=4202

% epsilon = 0.1|xMax – xMin|
epsilon = 0.1 * maxDiam;

recurrpt = tdrecurr_y(recurdata,epsilon, plotFigure);% plotFigure

% check it Marwan's method is equal to mine:
recur_plot_marwan = crp(sig,[],eDim,tau,epsilon,'euclidean','nonorm');
recur_plot = recurdata<=epsilon;
check = isequal(recur_plot_marwan,recur_plot);

% Marwan's trapping time:
tt_marwan = tt(recur_plot);
% it is similar to mine... the small difference might be the result of
% considering the recurr_plot_no_LOI instead of the whole recurr_plot 


%Recurrence quantification analysis
% rqa_stat - RQA statistics - [REC DET LMAX L ENT LAM TT]
% The main advantage of the recurrence quantification analysis is that it 
% can provide useful information even for short and non-stationary data (wiki)
% rqa_stat = recurrqa_y(recurrpt); % original
% the difference is in ENT, REC and the new measure L
rqa_stat = recurrence_quantification_analysis(recurrpt);
rqa_stat(end) = tt_marwan;

end