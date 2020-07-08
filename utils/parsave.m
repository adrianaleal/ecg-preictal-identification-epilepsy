function parsave(pathname, variables, variables_name, variable_name2)

if nargin<4
    variable_name2 = variables_name;
end



if numel(variables)>1
    string2eval = 'save(fullfile(pathname,[variable_name2 ''.mat''])';
    for ii = 1:numel(variables)
        eval([variables_name{ii} ' = variables{ii};'])
        string2eval = [string2eval ', ''' variables_name{ii} ''''];
    end
    string2eval = [string2eval ')'];
    eval(string2eval)
else
    eval([variables_name ' = variables{:};'])
    save(fullfile(pathname,[variable_name2 '.mat']), variables_name)
end
    
end
