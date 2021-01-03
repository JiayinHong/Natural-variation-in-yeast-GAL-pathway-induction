% define global variables here

perc_to_nm = 5.6*10^7;  % CONSTANT, the concentration of 1% glucose or galactose = 5.6*10^7 nM

% the GAL1 level under the highest gal/glu ratio corresponds fully ON state GAL1 level
logyfp_to_nm = @(x) exp(x-7.34) * 3000;
% exp(7.34) is the average fluorescent intensity while full induction