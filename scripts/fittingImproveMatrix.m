% this script displays the matrix that how much tuning each parameter
% one-at-a-time can improve the fitting of 15 natural yeast isolates
% created by JH on 2020.06.16, for Fig 4

% load obj results from each file, then store them all in 'obj_all'
clear
tmp = dir('../paramScanResults/');
tmp = tmp([tmp.isdir]~=1);  % filter out the directory
obj_all = struct;
for i=1:length(tmp)
    if ~strcmp(tmp(i).name, '.DS_Store')
        load(fullfile('../paramScanResults/',tmp(i).name))
        fdn = fieldnames(obj);
        fdn = fdn{1};
        obj_all.(fdn) = obj.(fdn);
    end
end

% store all the strain names into a cell array
clear strNames
nStr = numel(fieldnames(obj_all));
strNames = cell(nStr,1);

for i=1:length(allTraits)
    tmp = allTraits(i).name;
    strNames{i} = tmp(1:end-4);   % to remove '.mat'
end

% organize the data
nPerturb = size(obj_all.str01,1);
mydata = nan(nStr,nPerturb);
for iStr = 1:nStr
    fdn = num2str(sprintf('str%02d',iStr));
    for iPerturb = 1:nPerturb
        % the values are loss functions, thus take the min as the best fit
        mydata(iStr,iPerturb) = nanmin(obj_all.(fdn)(iPerturb,:));
    end
    % add a column to store the obj value without any perturbation
    % 24 perturbations flanking the original params
    mydata(iStr,nPerturb+1) = obj_all.(fdn)(1,25);
end
% compute how much loss reduced
diff_matrix = (mydata(:,end) - mydata(:,1:end-1)) ./ mydata(:,end);

load('../metadata/best_param.mat','param')
param_names = fieldnames(param);
rowLabels = param_names;
colLabels = strNames;

% make heatmap
CT = cbrewer('seq', 'Purples', 8);
figure
set(gcf,'position',[36 256 1231 442])
hm = heatmap(rowLabels, colLabels, diff_matrix, 'Colormap', CT, 'ColorLimits', [0 1]);
set(gca,'FontSize',16,'FontName','TimesNewRoman')
% export_fig('../figures/fitImproveMatrix','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')

% below show numbers of each element
% hm.CellLabelFormat = '%.1f';
% export_fig('../figures/fitImproveMatrix_show_data','-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')

%% print most improved fitting in terms of [strain,parameter]
thresh = 0.7;   % set threshold=0.7 to be improved significantly
[row,col] = find(hm.ColorDisplayData > thresh);
fileID = fopen('exp.txt','w');
fprintf(fileID,'[strain,params] that improved more than %.1f\n\n',thresh);
for i=1:numel(row)
    fprintf(fileID,'[%s\t%s\t%.2f]\n',colLabels{row(i)},param_names{col(i)},hm.ColorDisplayData(row(i),col(i)));
end
fclose(fileID);

