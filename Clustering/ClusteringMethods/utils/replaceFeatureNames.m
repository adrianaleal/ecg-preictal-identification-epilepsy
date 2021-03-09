function feat_names2analyse = replaceFeatureNames(feat_names2analyse)

% Replace the names of the features by the ones that should be plotted

feat_names2analyse = regexprep(feat_names2analyse,'_',' ');

index = find(~cellfun(@isempty,strfind(feat_names2analyse,'alpha')));
if ~isempty(index)
    for ll = 1:numel(index)
        % feat_names2analyse(index(ll)) = {[strrep(feat_names2analyse{index(ll)},' alpha',' $\alpha_') '$']}; % latex
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},' alpha',' \alpha_')};
    end

end

index = find(strcmp(feat_names2analyse,'SD1'));
if ~isempty(index)
    for ll = 1:numel(index)
        % feat_names2analyse(index(ll)) = {[strrep(feat_names2analyse{index(ll)},'SD','SD$_') '$']}; % latex
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'SD','SD_')};
    end
end

index = find(strcmp(feat_names2analyse,'SD2'));
if ~isempty(index)
    for ll = 1:numel(index)
        % feat_names2analyse(index(ll)) = {[strrep(feat_names2analyse{index(ll)},'SD','SD$_') '$']}; % latex
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'SD','SD_')};
    end
end

index = find(strcmp(feat_names2analyse,'SD1toSD2'));
if ~isempty(index)
    for ll = 1:numel(index)
        % feat_names2analyse(index(ll)) = {[strrep(feat_names2analyse{index(ll)},'SD1toSD2','SD$_1$/SD$_2') '$']}; % latex
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'SD1toSD2','SD_1/SD_2')};
    end
end


index = find(strcmp(feat_names2analyse,'LFtoHF'));
if ~isempty(index)
    for ll = 1:numel(index)
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'LFtoHF','LF/HF')};
    end
end


index = find(strcmp(feat_names2analyse,'largLyapunov'));
if ~isempty(index)
    for ll = 1:numel(index)
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'largLyapunov','LLE')};
    end
end

index = find(strcmp(feat_names2analyse,'corrDim'));
if ~isempty(index)
    for ll = 1:numel(index)
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'corrDim','CD')};
    end
end

index = find(strcmp(feat_names2analyse,'RQA_Lmax'));
if ~isempty(index)
    for ll = 1:numel(index)
        feat_names2analyse(index(ll)) = {strrep(feat_names2analyse{index(ll)},'RQA_Lmax','RQA L_{max}')};
    end
end

end


