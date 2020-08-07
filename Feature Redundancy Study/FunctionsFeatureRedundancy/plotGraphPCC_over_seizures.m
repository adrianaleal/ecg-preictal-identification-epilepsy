function [h] = plotGraphPCC_over_seizures(x, G)


g_nodes = table2cell(G.Nodes);

% pentagon
L = linspace(0,2*pi,6);
xv = 1.5*cos(L)';
yv = 1.5*sin(L)';
% RRVar, TOTAL POWER, SD2, SDNN, VLF POWER

variables = {'RRMax', 'RRMean', 'RRmin', 'RQA L', 'RQA TT', 'RQA ENT', ...
    'RQA DET', 'RQA LAM', 'SampEn', 'ApEn', 'DFA alpha2', 'LF NORM', ...
    'HF NORM', 'SD1toSD2', 'pNN50', 'SD1', 'RMSSD', 'SDSD', 'NN50', ...
    'HF POWER', 'RRVar', 'TOTAL POWER', 'SD2', 'SDNN', 'VLF POWER'};
positionsx = [8 8 8 5 6 5.5 5.5 5.5 7 7.5 7 6 6.5 6.5 2 2 1 1 3 ...
    3 xv(1:5)'+2];

positionsy = [2 3 1 0 0 1 2 3 5 5 9 9 8 7 0 4 1 3 1 3 yv(1:5)'+8];

x_changed = zeros(size(x));
y_changed = x_changed;
for nn = 1:numel(g_nodes)
    ind_variable = strcmp(g_nodes(nn), variables);
    x_changed(nn) = positionsx(ind_variable);
    y_changed(nn) = positionsy(ind_variable);
end

h = plot(G,'XData', x_changed, 'YData', y_changed);

feat_names = replaceFeatureNames(h.NodeLabel);


down_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'pNN50', 'VLF POWER', 'SD1toSD2'}); 
down_feat(indexes) = 1;
down_feat = logical(down_feat);
text(h.XData(down_feat), h.YData(down_feat)-0.2 ,feat_names(down_feat), ...
    'VerticalAlignment','top',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)


indexes = ismember(g_nodes, {'TOTAL POWER', 'RQA LAM', 'SD1', 'SampEn'});
up_feat = zeros(numel(x),1);
up_feat([indexes]) = 1;
up_feat = logical(up_feat);
text(h.XData(up_feat), h.YData(up_feat)+0.2,feat_names(up_feat), ...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)


left_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'RMSSD', 'SDSD', 'SD2', 'SDNN', 'RQA DET', 'RQA ENT'});
left_feat(indexes) = 1;
left_feat = logical(left_feat);
text(h.XData(left_feat)-0.2, h.YData(left_feat),feat_names(left_feat), ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'right',...
    'FontSize', 10)

down_left_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'RQA L'});
down_left_feat(indexes) = 1;
down_left_feat = logical(down_left_feat);
text(h.XData(down_left_feat)-0.2, h.YData(down_left_feat)-0.2, ...
    feat_names(down_left_feat), ...
    'VerticalAlignment','top',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)

down_right_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'RQA TT'});
down_right_feat(indexes) = 1;
down_right_feat = logical(down_right_feat);
text(h.XData(down_right_feat)+0.2, h.YData(down_right_feat)-0.2, ...
    feat_names(down_right_feat), ...
    'VerticalAlignment','top',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)


up_left_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'LF NORM'});
up_left_feat(indexes) = 1;
up_left_feat = logical(up_left_feat);
text(h.XData(up_left_feat)-0.2, h.YData(up_left_feat)+0.2, ...
    feat_names(up_left_feat), ...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)


up_right_feat = zeros(numel(x),1);
indexes = ismember(g_nodes, {'DFA alpha2'});
up_right_feat(indexes) = 1;
up_right_feat = logical(up_right_feat);
text(h.XData(up_right_feat)+0.2, h.YData(up_right_feat)+0.2, ...
    feat_names(up_right_feat), ...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment', 'center',...
    'FontSize', 10)

other_feat = ~sum([up_feat, left_feat, down_feat, down_left_feat, ...
    up_left_feat, down_right_feat, up_right_feat],2);
text(h.XData(other_feat)+0.2, h.YData(other_feat), feat_names(other_feat), ...
    'VerticalAlignment','Middle',...
    'HorizontalAlignment', 'left',...
    'FontSize', 10)

end