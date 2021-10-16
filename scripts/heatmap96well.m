% this script displays comparison between simulation results
% and experimental data
% created by JH on 2020.06.08, for Fig 2 & Fig 5
% updated by JH on 2021.05.31, add GAL3 swap experiment

%% for Fig 2, showing the fitting for S288C
for dataType = {'wildtype', 'mig1d', 'gal80d'}
    dataType = dataType{1};
    switch dataType
        case 'wildtype'
            load('../metadata/best_param.mat')
            plot_heatmap(param, dataType)
        case 'mig1d'
            load('../metadata/best_param.mat')
            param.aR = 0;
            plot_heatmap(param, dataType)
        case 'gal80d'
            load('../metadata/best_param.mat')
            param.a80 = 0;
            param.ag80 = 0;
            plot_heatmap(param, dataType)
    end
end

%% for Fig 5, showing phenotypic switch
for strains = {'I14','L-1528','273614N','YPS163'...
        ,'YPS606','IL-01','UWOPS87','YPS606_GAL3_swap'}
    strain = strains{1};
    load('../metadata/best_param.mat')
    switch strain
        case 'I14'
            param.rHXT = param.rHXT * 0.0075;
        case 'L-1528'
            param.kf83 = param.kf83 * 75;
        case '273614N'
            param.kf83 = param.kf83 * 100;
        case 'YPS163'
            param.rHXT = param.rHXT * 0.013; 
        case 'YPS606'
            param.kf3 = param.kf3 * 42;
        case 'IL-01'
            param.kf83 = param.kf83 * 10;
        case 'UWOPS87'
            param.rHXT = param.rHXT * 0.02;
        case 'YPS606_GAL3_swap'
    end 
    plot_heatmap(param,strain)
end

function plot_heatmap(param, dataType)
%% claim saving directory
saveDir = fullfile(sprintf('../96wellPlot/%s', dataType), datestr(now,'yy-mm-dd'));  % the directory to store the figures
if ~isdir(saveDir)
    mkdir(saveDir)
end
%% load experimental data
switch dataType
    case 'wildtype'
        load('../traits/13-S288C.mat')
    case 'mig1d'
        load('../metadata/mig1d_all_data.mat')
    case 'gal80d'
        load('../metadata/gal80d_all_data.mat')
    case 'DBVPG1853'
        load('../traits/01-DBVPG1853.mat')
    case 'I14'
        load('../traits/02-I14.mat')  
    case 'IL-01'
        load('../traits/03-IL-01.mat')
    case 'L-1528'
        load('../traits/04-L-1528.mat')
    case 'YPS163'
        load('../traits/05-YPS163.mat')
    case 'YPS606'
        load('../traits/06-YPS606.mat')
    case '273614N'
        load('../traits/07-273614N.mat')
    case 'BC187'
        load('../traits/08-BC187.mat')
    case 'CLIB324'
        load('../traits/09-CLIB324.mat')
    case 'DBVPG6765'
        load('../traits/10-DBVPG6765.mat')
    case 'UWOPS87'
        load('../traits/11-UWOPS87-2421.mat')
    case 'Y12'
        load('../traits/12-Y12-SGRP.mat')
    case 'YJM975'
        load('../traits/14-YJM975.mat')
    case 'YJM981'
        load('../traits/15-YJM981.mat')
    case 'YPS606_GAL3_swap'
        load('../traits/16-GAL3swap.mat')
    otherwise
        error('strain records not found')
end

expt_96well = trait;
galLabel = {'None','2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1','2','4'};
gluLabel = {'None','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1'};
load_global
alldata = nan(8,12);

if strcmp(dataType, 'gal80d')
    % gal80d is different from the other two, because all the mask_basal in
    % gal80d equals to 0 (unimodal), thus there's no need to replace
    % induced level using basal level, just change the wells we don't trust
    % to NaN
    ind3 = expt_96well.mask_induction == 0;
    expt_96well{ind3, 'ind_level'} = NaN;
else
    % to clean the data
    ind1 = find(expt_96well.mask_induction == 0);     % all the rows whose mask_induction == 0
    tmp = expt_96well(ind1,:).mask_basal == 0;
    ind1(tmp) = [];
    
    % use basal_level to represent induced level in ind1
    expt_96well(ind1,:).ind_level = expt_96well(ind1,:).basal_level;
end

alldata(:) = logyfp_to_nm(expt_96well{:,'ind_level'});

%% fetch simulation results
output = evalGalPathway(param, trait, '96well');

simG1_96well = nan(8,12);
simG1_96well(:) = output.all_conc_Gal(:,1);

if exist('ind1', 'var')
    % use basal level to compare with expt data when mask_basal = 1 &&
    % mask_induction == 0
    simG1_96well(ind1) = output.all_conc_Glu(ind1,1);
end

tmp_exp = log(alldata);
max_exp = max(tmp_exp(:));
tmp_sim = log(simG1_96well);

% first, heatmap for the expt trait
figure
if ~strcmp(dataType,'YPS606_GAL3_swap')
    set(gcf, 'position', [281 117 383 581])
    subplot(2,1,1)
else
    set(gcf,'position',[281 389 383 309])
end
% to cope with missing values in experimental data
% i.e. leave the color as white instead of using the lowest value to pad
imAlpha = ones(size(tmp_exp));
imAlpha(isnan(tmp_exp)) = 0;
h1 = imagesc(tmp_exp./max_exp,'AlphaData',imAlpha);
% only missing values in expt. not in sims
h1.Parent.XTick = 1:12;
h1.Parent.YTick = 1:8;
h1.Parent.XTickLabel = galLabel;
h1.Parent.YTickLabel = fliplr(gluLabel);
h2=colorbar;
if strcmp(dataType, 'gal80d')
    h2.Limits = [0.6 1];
    h2.Ticks = [0.6, 0.7, 0.8, 0.9, 1.0];
else
    h2.Limits = [0 1];
    h2.Ticks = [0, 0.2, 0.4, 0.6, 0.8, 1];
end
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',20,'FontName','Times New Roman','FontWeight','normal')
colorLim = get(gca,'CLim');
tt = title(sprintf('%s expt G1 induction level', dataType));

% second, heatmap for simulation results
if ~strcmp(dataType,'YPS606_GAL3_swap')
    subplot(2,1,2)
    h3=imagesc(tmp_sim./max_exp);
    h3.Parent.XTick = 1:12;
    h3.Parent.YTick = 1:8;
    h3.Parent.XTickLabel = galLabel;
    h3.Parent.YTickLabel = fliplr(gluLabel);
    h4=colorbar;
    h4.Limits = h2.Limits;  % set the same color limits for colorbar
    h4.Ticks = h2.Ticks;
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',20,'FontName','Times New Roman','FontWeight','normal')
    set(gca,'CLim',colorLim)  % set the same color limits for heatmap wells
    tt = title(sprintf('%s sims G1 induction level', dataType));
    % export_fig(fullfile(saveDir, sprintf('%s induction level', dataType)), '-pdf', '-transparent', '-c[NaN,NaN,NaN,NaN]')

end
end