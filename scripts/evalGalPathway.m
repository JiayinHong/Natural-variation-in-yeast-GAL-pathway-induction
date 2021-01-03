function output = evalGalPathway(param, trait, fit_type)
%   the third version of eval_param, aim to increase the compatibility
%   when dealing with multiple rows and columns of data
fit_type_config;
trait = trait(index_list,:);

output = getInitParamsGalPathway(param);
y_ss_Glu = evalMultiSugarConcentrations( param, output.y0_Glu, trait.gluc, trait.galc );
y_ss_Gal = evalMultiSugarConcentrations( param, output.y0_Gal, trait.gluc(end:-1:1), trait.galc(end:-1:1) );
y_ss_Gal = y_ss_Gal(end:-1:1,:);

basal_level = y_ss_Glu(:,1);
induction_level = y_ss_Gal(:,1);

[output.G1obj, ~, ~] = calculate_obj( trait, basal_level, induction_level );
output.all_conc_Glu = y_ss_Glu; % all 12 variables concentration at steady state, initial from Glu only condition
output.all_conc_Gal = y_ss_Gal; % all 12 variables concentration at steady state, initial from Gal only condition

% for GAL1 level, add autofluorescence
autofluorescence = get_auto_fluorescence(trait);
output.all_conc_Glu(:,1) = output.all_conc_Glu(:,1) + autofluorescence;
output.all_conc_Gal(:,1) = output.all_conc_Gal(:,1) + autofluorescence;

end
