function [ sum_obj, simulation_result_linear, experiment_result_linear ] = calculate_obj( trait, basal_level, induction_level )
% this function takes in the basal level and induction level of the simulation results,
% and the experimental data (trait), then calculate the difference between them,
% return the sum of the difference as sum_obj (the objective function)

load_global;
autofluorescence = get_auto_fluorescence(trait);

expt_basal = trait.basal_level;
expt_induction = trait.ind_level;

% check if there is nan in basal or induced
tmp = isnan(expt_basal);
expt_basal(tmp) = expt_induction(tmp);
tmp = isnan(expt_induction);
expt_induction(tmp) = expt_basal(tmp);

% mask, expt and simulation all have two columns,
% the first - basal, and the second - induced
mask = [trait.mask_basal, trait.mask_induction];
experiment_result = [expt_basal, expt_induction];
simulation_result = [basal_level, induction_level];

simulation_result_linear = autofluorescence + simulation_result;
% convert the expt result from natural log to linear
experiment_result_linear = logyfp_to_nm(experiment_result);

simulation_result = log10(simulation_result_linear);
experiment_result = log10(experiment_result_linear);

sum_obj = (simulation_result - experiment_result).^2 .* mask;

sum_obj = nansum(sum_obj(:));

end
