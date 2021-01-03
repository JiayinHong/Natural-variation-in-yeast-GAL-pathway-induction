function [ autofluorescence ] = get_auto_fluorescence( trait )
%   define the minimum of the basal level as auto-fluorescence
load_global
autofluorescence = logyfp_to_nm(min(trait.basal_level));

end