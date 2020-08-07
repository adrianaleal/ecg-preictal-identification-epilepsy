function pats_name = getPatientsNames(patients_folder_path)

patients = dir(patients_folder_path);
pats_name = {patients.name}';
split1 = regexp(pats_name, '_', 'split');
get_pat_ids = vertcat(split1{:});
[~,I] = sort(cellfun(@str2double, get_pat_ids(:,2)));
pats_name = pats_name(I);

end