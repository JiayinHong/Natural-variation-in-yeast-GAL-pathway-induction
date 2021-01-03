% this script clustering parameter perturbations based on their effects on
% bimodality and expression level, then use t-SNE to visualize parameter
% clusters
% created by JH on 2020.06.11, for Fig 3 & Supplementary Fig 2-5

% store the start and end concentration of galactose that exhibit bimodality,
% as well as the max level for a given gluc conc.
% every 3 cols are of the same gluc conc
load('../metadata/best_param.mat')
load('../traits/13-S288C.mat')
output = evalGalPathway(param, trait, '96well');
default_on = output.all_conc_Gal;
default_off = output.all_conc_Glu;

galTitrate = 1:8:89;
galGrad = {'None','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2'};
gluGrad = {'0','-1','-2','-3','-4','-5','-6','None'};

default_stat = [];
for i_row = 1:8
    index = i_row-1+galTitrate;
    tmp = [default_on(index) > 5* default_off(index)]';
    row_vec(1) = find(tmp==1,1,'first');    % 1 col is the start well for bimodal
    row_vec(2) = find(tmp==1,1,'last');     % 2 col is the end well for bimodal
    row_vec(3) = max(default_on(index));    % 3 col is the max level for the fixed gluc
    default_stat = [default_stat,row_vec];
end

%% calculate and store the features of a given parameter with a specific perturbation
load('../metadata/param_perturb_Fig3.mat')

nParam = size(ONpeak,1);
nParam = nParam - 9;    % to exempt Hill coefficients
nPerturb = size(ONpeak,2);

all_matrix = [];
for iParam = 1:nParam
    for iPerturb = 1:nPerturb
        row_vec_all = [];
        for i_row = 1:8
            index = i_row-1+galTitrate;
            cur_ON = ONpeak{iParam,iPerturb};
            cur_OFF = OFFpeak{iParam,iPerturb};
                        
            row_vec(3) = max(cur_ON(index));    % 3 col is the max level for the fixed gluc
            
            tmp = [cur_ON(index) > 5* cur_OFF(index)]';
            if tmp==0   % unimodal
                if row_vec(3) > max(ONpeak{iParam,iPerturb}(index))
                    row_vec(1:2) = 0;   % always ON
                else
                    row_vec(1:2) = 13;  % always OFF
                end
            else
                row_vec(1) = find(tmp==1,1,'first');    % 1 col is the start well for bimodal
                row_vec(2) = find(tmp==1,1,'last');     % 2 col is the end well for bimodal
            end
            row_vec_all = [row_vec_all,row_vec];
        end
        
        row_vec_all(1:3:22) = row_vec_all(1:3:22) - default_stat(1:3:22);
        row_vec_all(2:3:23) = row_vec_all(2:3:23) - default_stat(2:3:23);
        row_vec_all(3:3:end) = (row_vec_all(3:3:end) - default_stat(3:3:end)) ./ default_stat(3:3:end);

        all_matrix = [all_matrix;row_vec_all];
    end
end

%% generate row labels and col labels
% to generate row labels (parameter name * perturbation factor)
param_names = fieldnames(param);
observations = cell(1,216);     % row labels
ind = 1;
for iParam = 1:nParam
    for perturb_factor = {'*0.01','*0.1','*0.5','*2','*10','*100'}
        tmp = perturb_factor{1};
        observations{ind} = [param_names{iParam},tmp];
        ind = ind+1;
    end
end
% to generate col labels (gluc conc + induction ON / full / max level)
ind = 1;
features = cell(1,24);      % col labels
for gluc = {'1%','0.5%','0.25%','0.125%','0.0625%','0.0312%','0.0156%','None'}
    gluc_conc = gluc{1};
    for attribute = {' ON',' full',' level'}
        tmp = attribute{1};
        if strcmp(tmp,' full')
            features{ind} = ['gluc=',gluc_conc,':',tmp];
        else
            features{ind} = tmp;
        end
        ind = ind+1;
    end
end

%% plot the cluster dendrogram
CG = clustergram(all_matrix, 'RowLabels',observations, 'ColumnLabels',features ...
        , 'Standardize','none', 'Cluster','column', 'Colormap',redbluecmap ...
        , 'ShowDendrogram','on', 'LogTrans', 0);
    
rm = struct('GroupNumber',{187,192,202},'Annotation',{'more inducible','bimodailty insensitive','more repressed'},...
     'Color',{'r',[1 1 0],'b'});
set(CG,'RowGroupMarker',rm)

cgAxes = plot(CG);
cgAxes.Position(2) = cgAxes.Position(2) + 0.15;
% Children(3) are CMarkerRowAxes
cgAxes.Parent.Children(3).Position(2) = cgAxes.Parent.Children(3).Position(2) + 0.15;
colorbar(cgAxes,'location','eastoutside')
set(gcf, 'position', [322 50 863 655])
set(cgAxes, 'XTickLabelRotation',45)
set(cgAxes, 'FontSize', 17, 'FontName', 'TimesNewRoman')
% export_fig('../figures/clustergram of perturbations regarding bimodality variation', '-pdf', '-transparent', '-c[NaN,NaN,NaN,NaN]')

%% plot group of parameters that make GAL more inducible
load('../metadata/group_inducible.mat')
CG187 = clustergram(group187.ExprValues, 'RowLabels',group187.RowNodeNames ...
    ,'ColumnLabels',group187.ColumnNodeNames, 'Standardize','none' ...
    ,'Cluster','column','Colormap',redbluecmap, 'LogTrans', 0, 'ShowDendrogram','off');
cgAxes = plot(CG187);
cgAxes.Position(2) = cgAxes.Position(2) + 0.1;
colorbar(cgAxes,'location','eastoutside')
set(cgAxes, 'FontSize', 14, 'FontWeight', 'normal')
set(cgAxes, 'XTickLabelRotation',45)
tithand = addTitle(CG187,'group ''more inducible''');
tithand.Position(2) = tithand.Position(2) - 0.05;
tithand.FontWeight = 'bold';
set(gcf,'position',[360 108 738 590])
% export_fig('../figures/group more inducible','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')

%% plot group of parameters that are insensitive to bimodality
load('../metadata/group_insensitive.mat')
CG192 = clustergram(group192.ExprValues, 'RowLabels',group192.RowNodeNames ...
    ,'ColumnLabels',group192.ColumnNodeNames, 'Standardize','none' ...
    ,'Cluster','column','Colormap',redbluecmap, 'LogTrans', 0, 'ShowDendrogram','off');
cgAxes = plot(CG192);
cgAxes.Position(2) = cgAxes.Position(2) + 0.1;
colorbar(cgAxes,'location','eastoutside')
set(cgAxes, 'FontSize', 14, 'FontWeight', 'normal')
set(cgAxes, 'XTickLabelRotation',45)
tithand = addTitle(CG192,'group ''bimodality insensitive''');
tithand.Position(2) = tithand.Position(2) - 0.05;
tithand.FontWeight = 'bold';
set(gcf,'position',[360 108 738 590])
% export_fig('../figures/group bimodality insensitive','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')

%% plot group of parameters that make GAL more repressed
load('../metadata/group_repressed.mat')
CG202 = clustergram(group202.ExprValues, 'RowLabels',group202.RowNodeNames ...
    ,'ColumnLabels',group202.ColumnNodeNames, 'Standardize','none' ...
    ,'Cluster','column','Colormap',redbluecmap, 'LogTrans', 0, 'ShowDendrogram','off');
cgAxes = plot(CG202);
cgAxes.Position(2) = cgAxes.Position(2) + 0.1;
colorbar(cgAxes,'location','eastoutside')
set(cgAxes, 'FontSize', 14, 'FontWeight', 'normal')
set(cgAxes, 'XTickLabelRotation',45)
tithand = addTitle(CG202,'group ''more repressed''');
tithand.Position(2) = tithand.Position(2) - 0.05;
tithand.FontWeight = 'bold';
set(gcf,'position',[360 108 738 590])
% export_fig('../figures/group more repressed','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')

%% tSNE visualization
group_more_inducible = group187.RowNodeNames;
group_bimodal_insensitive = group192.RowNodeNames;
group_more_repressed = group202.RowNodeNames;

% generate a cell array to store the group information
groups = cell(216,1);
groups(ismember(observations, group_more_inducible)) = {'more inducible'};
groups(ismember(observations, group_bimodal_insensitive)) = {'bimodality insensitive'};
groups(ismember(observations, group_more_repressed)) = {'more repressed'};
for i=1:numel(groups)
    if isempty(groups{i})
        groups(i) = {'others'};
    end
end

figure;
set(gcf,'position',[360 89 836 609]);
[Y,loss] = tsne(all_matrix,'Perplexity',50,'Algorithm','exact','Distance','Euclidean');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),groups)
title(sprintf('Euclidean loss=%.3f',loss))
set(gca,'FontSize',14)

[Y,loss] = tsne(all_matrix,'Perplexity',50,'Algorithm','exact','Distance','Chebychev');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),groups)
title(sprintf('Chebychev loss=%.3f',loss))
set(gca,'FontSize',14)

[Y,loss] = tsne(all_matrix,'Perplexity',50,'Algorithm','exact','Distance','Minkowski');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),groups)
title(sprintf('Minkowski loss=%.3f',loss))
set(gca,'FontSize',14)

[Y,loss] = tsne(all_matrix,'Perplexity',50,'Algorithm','exact','Distance','Cityblock');
subplot(2,2,4)
gscatter(Y(:,1),Y(:,2),groups)
title(sprintf('Cityblock loss=%.3f',loss))
set(gca,'FontSize',14)

% export_fig('../figures/tSNE visualization','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')


