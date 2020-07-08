function [variable] = parload(pathname, variable_name_out, ...
    variable_name_in)


if nargin==2
    variable_name_in = variable_name_out;
end

load(fullfile(pathname,[variable_name_out '.mat']), variable_name_in)
eval(['variable = ' variable_name_in ';'])

end