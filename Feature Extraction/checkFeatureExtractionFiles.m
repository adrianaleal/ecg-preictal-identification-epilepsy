function check = checkFeatureExtractionFiles(pat_name, patients_info_final, ...
    features_folder_path)


sub_clustering_folder_path = 'ExtractedFeaturesNoEctopicBeatsOverlap';

foldersInPath = dir(fullfile(features_folder_path, sub_clustering_folder_path));
foldersInPath = foldersInPath(~ismember({foldersInPath.name},{'.','..'}));

folders_name = {foldersInPath.name}';
split1 = regexp(folders_name, '_', 'split');

get_folder_ids = vertcat(split1{:});
[~,I] = sort(cellfun(@str2double,get_folder_ids(:,2)));
folder_name = folders_name(I);


save_info = cell(numel(pat_name),1);
save_info_date = save_info;
results_info = save_info;
results_info_date = save_info;
save_info_intervals = save_info;
save_info_intervals_date = save_info;

for pp = 1:numel(pat_name)
    
    folderPath = fullfile(features_folder_path, sub_clustering_folder_path, pat_name{pp});
    
    if exist(folderPath, 'dir')
        filesInFolder = dir(folderPath);
        
        seiz_folders = filesInFolder(~ismember({filesInFolder.name},{'.','..'}));
        n_seiz = numel(seiz_folders);
        n_seiz2 = patients_info_final(pp).n_seizures- ...
            patients_info_final(pp).n_discarded_seizures_by_no_file- ...
            patients_info_final(pp).n_discarded_seizures_by_240min;
        
        seiz_save_info = zeros(n_seiz2,1);
        
        seiz_save_info_date = cell(n_seiz2,1);
        
        for ss = 1:n_seiz
            seizFolderPath = fullfile(folderPath,seiz_folders(ss).name);
            filesInSeizureFolder = dir(seizFolderPath);
            filesInSeizureFolder = filesInSeizureFolder(~ismember({filesInSeizureFolder.name},{'.','..'}));
            
            name = [pat_name{pp} '_seizure' num2str(ss)];
            
            ind = strcmp(vertcat({filesInSeizureFolder.name}), [name '.mat']);
            if any(ind)
                seiz_save_info_date(ss) = {filesInSeizureFolder(ind).date};
            end
            if datetime(filesInSeizureFolder(ind).date, 'Format','dd-MMM-yyyy HH:mm:ss')> ...
                    datetime('14-Jan-2020 09:00:00', 'Format','dd-MMM-yyyy HH:mm:ss')
                seiz_save_info(ss) = 1;
            end
        end
        
        if all(seiz_save_info==1)
            save_info(pp) = {1};
        elseif all(seiz_save_info==0)
            save_info(pp) = {0};
        else
            save_info(pp) = {seiz_save_info};
        end
        
        
        if any(seiz_save_info)
            save_info_date(pp) = {seiz_save_info_date};
        end
        
    end
    
end

check = [pat_name, save_info, save_info_date];

end
