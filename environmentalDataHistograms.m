% This script reads in the 'plantFoodsData.mat' dataset and the 'bird_withinChannelTbl.mat' table, the 
% latter of which is made by the 'build_tables.m script. From these, it extracts data on sun elevation, 
% sun occlusion, and cloud cover for each image and object comparison and plots their distributions in 
% histograms. 

%% PER IMAGE
clear
close all
% load data
load('plantFoodsData.mat');

% preallocation
numImgSets = length(data);
sunElev = zeros(numImgSets, 1);
sunOccl = zeros(numImgSets, 1); 
cloudCover = zeros(numImgSets, 1); 

% extract environmental data from full dataset
for i = 1:numImgSets
    sunElev(i) = str2double(data(i).sun_elevation);
    if isempty(data(i).sun_occlusion)
        sunOccl(i) = nan;
    else
        sunOccl(i) = data(i).sun_occlusion;
    end
    if isempty(data(i).cloud_cover)
        cloudCover(i) = nan;
    else
        cloudCover(i) = data(i).cloud_cover;
    end
end

% plot histograms
histdata = table(sunElev, sunOccl, cloudCover);
plotHistograms(histdata, 'Per image', '# of images')

%% PER OBJECT COMPARISON
load('bird_withinChannelTbl.mat');
Tbl = withinChTable;
clear withinChTable

% limits object comparisons to objects under similar illumination
Tbl = sameIllum(Tbl);

% reclassifies 'grass upper' as 'upper leaf surfaces' and 'grass lower' as 'lower leaf surfaces'
Tbl = reclassGrass(Tbl);

% reclassify fiddleheads as leaf buds
Tbl = reclassFiddleheads(Tbl);

% reclassify 'unripe fleshy fruit capsule' as 'unripe fleshy fruit' 
Tbl = reclassUnripeFFCapsules(Tbl);

%reclassify various botanical names for immature non-fleshy fruits as 'immature seed'
Tbl = reclassImmatSeeds(Tbl);

% object categories to compare
xObjCats = {'upper leaf surfaces';
               'lower leaf surfaces';
               'woody branch or trunk'};
yObjCats = {'leaf buds';
            'new leaves - upper';
            'new leaves - lower';
            'immature seed';
            'unripe fleshy fruit'};

% organizes all object comparisons into tables and plots histograms
fullTbl = table();
for c = 1:size(yObjCats,1)
    for i = 1:size(xObjCats,1) 
        rows = (Tbl.obj1_class == xObjCats{i} & Tbl.obj2_class == yObjCats{c}) | ...
               (Tbl.obj2_class == xObjCats{i} & Tbl.obj1_class == yObjCats{c}); 
        selectTbl = Tbl(rows,:);
        if height(selectTbl) > 5
            fullTbl = vertcat(fullTbl, selectTbl);
            plotHistograms(selectTbl, ['Per ', yObjCats{c}, ' vs ', xObjCats{i}], '# obj. comps.')
        end
    end
end
plotHistograms(fullTbl, 'Per all object comparisons', '# obj. comps.')

%% FUNCTIONS
% plots data
function [] = plotHistograms(fullTbl, food, yAxLabel)
    % sun elevation
    figure
    subplot(1,3,1)
    histogram(fullTbl.sunElev, 'BinEdges', (-5:5:80), 'FaceColor', 'k')
    set(gca,'TickDir','out')
    set(gca, 'FontSize', 9)
    ylabel(yAxLabel)
    xlabel('sun elevation')

    % sun occlusion
    subplot(1,3,2)
    histogram(fullTbl.sunOccl, 'BinEdges', (0.5:5.5), 'FaceColor', 'k')
    set(gca,'TickDir','out')
    xticklabels([0 25 50 75 100])
    set(gca, 'FontSize', 9)
    ylabel(yAxLabel)
    xlabel('sun occlusion')
    title(food)

    % cloud cover
    subplot(1,3,3)
    histogram(fullTbl.cloudCover, 'BinEdges', (0.5:5.5), 'FaceColor', 'k')
    set(gca,'TickDir','out')
    xticklabels([0 25 50 75 100])
    set(gca, 'FontSize', 9)
    set(gcf, 'Position', [1 1 420 90])
    ylabel(yAxLabel)
    xlabel('cloud cover')
end

% limits object comparisons to objects under similar illumination
function Tbl = sameIllum(Tbl)
    rows = (Tbl.directLightObj1 == '0' & Tbl.directLightObj2 == '0') | (Tbl.directLightObj1 == '1' & Tbl.directLightObj2 == '1') | (isundefined(Tbl.directLightObj1) & isundefined(Tbl.directLightObj2));
    Tbl = Tbl(rows,:);
end

% reclassifies 'grass upper' as 'upper leaf surfaces' and 'grass lower' as 'lower leaf surfaces'
function Tbl = reclassGrass(Tbl)
    rows = Tbl.obj1_class == 'grass upper';
    Tbl.obj1_class(rows) = 'upper leaf surfaces';
    rows = Tbl.obj2_class == 'grass upper';
    Tbl.obj2_class(rows) = 'upper leaf surfaces';
    
    rows = Tbl.obj1_class == 'grass lower';
    Tbl.obj1_class(rows) = 'lower leaf surfaces';
    rows = Tbl.obj2_class == 'grass lower';
    Tbl.obj2_class(rows) = 'lower leaf surfaces';
end

% reclassify fiddleheads as leaf buds
function Tbl = reclassFiddleheads(Tbl)
    rows = Tbl.obj1_class == 'fiddlehead';
    Tbl.obj1_class(rows) = 'leaf buds';

    rows = Tbl.obj2_class == 'fiddlehead';
    Tbl.obj2_class(rows) = 'leaf buds';
end

% reclassify unripe fleshy fruit capsules as unripe fleshy fruits
function Tbl = reclassUnripeFFCapsules(Tbl)
    rows = Tbl.obj1_class == 'unripe fleshy fruit capsule';
    Tbl.obj1_class(rows) = 'unripe fleshy fruit';

    rows = Tbl.obj2_class == 'unripe fleshy fruit capsule';
    Tbl.obj2_class(rows) = 'unripe fleshy fruit';
end

% reclassify various botanical names for immature non-fleshy fruits as 'immature seed'
function Tbl = reclassImmatSeeds(Tbl)
    rows = Tbl.obj1_class == 'immature achene' | Tbl.obj1_class == 'immature nut' | Tbl.obj1_class == 'immature seedpod' | Tbl.obj1_class == 'immature spikelet';
    Tbl.obj1_class(rows) = 'immature seed';

    rows = Tbl.obj2_class == 'immature achene' | Tbl.obj2_class == 'immature nut' | Tbl.obj2_class == 'immature seedpod' | Tbl.obj2_class == 'immature spikelet';
    Tbl.obj2_class(rows) = 'immature seed';
end