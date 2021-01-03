% this script is used to display the phenotypic diversity in
% natural yeast isolates
% updated by JH on 2020.06.04, for Fig 1

galLabel = {'None','2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1','2','4'};
gluLabel = {'None','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1'};
load_global
Cmap = parula;
strains = dir('../traits/');
strains = strains([strains.isdir]~=1,:);
if strcmp(strains(1).name,'.DS_Store')
    strains(1) = [];
end
emptyCell1 = cell(1,12);
emptyCell1(:) = {''};
emptyCell2 = cell(1,8);
emptyCell2(:) = {''};

figure
ha = tight_subplot(4,6,[.05 .03],[.05 .05],[.05 .05]);

set(gcf,'position',[162 77 1014 621])
for i=1:length(strains)
    axes(ha(i))
    if i==14 || i==15
        axes(ha(i+5))   % manually shift the row
    end
    strainName = strains(i).name;
    strainName = strainName(4:end-4);
    
    load(fullfile('../traits/',strains(i).name),'trait')
    expt_96well = trait;
    alldata = nan(8,12);
    
    % to clean the data
    ind1 = find(expt_96well.mask_induction == 0);     % all the rows whose mask_induction == 0
    tmp = find(expt_96well(ind1,:).mask_basal == 0);
    ind2 = ind1(tmp);         % the rows whose mask_basal also equals to 0
    ind1(tmp) = [];           % remove ind2 from ind1, so that ind1 only contains
                              % rows whose mask_induction == 0 while mask_basal ~= 0
    
    % use NaN for those induced level in ind2
    expt_96well{ind2, 'ind_level'} = NaN;
    % use basal_level to represent induced level in ind1
    expt_96well(ind1,:).ind_level = expt_96well(ind1,:).basal_level;
    
    alldata(:) = logyfp_to_nm(expt_96well{:,'ind_level'});
    tmp = log(flip(alldata));
    tmp1 = tmp./max(tmp(:));
    hp = pcolor(padarray(tmp1,[1 1],nan,'post'));
%     set(hp,'edgecolor','none')  % whether show grid or not

    if i==6
        hc=colorbar;    % then manually adjust the position of colorbar
        hc.Position(1) = hc.Position(1) + 0.04;
        hc.FontSize = 12;
        hc.FontName = 'TimesNewRoman';
    end
    set(gca,'XTickLabel','','YTickLabel','')
    
    if i==14
        set(gca,'XTick',[3:2:13],'YTick',1:1:8)
        set(gca,'XTickLabel',galLabel([2:2:12]),'YTickLabel',gluLabel(1:1:8))
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize', 12, 'FontName','TimesNewRoman')
    end
    axis('square')
    title(sprintf('%s', strainName),'FontSize',15,'FontName','TimesNewRoman');
end
hold all

% export_fig('../figures/natural_variation', '-pdf', '-transparent', '-c[NaN,NaN,NaN,NaN]')
