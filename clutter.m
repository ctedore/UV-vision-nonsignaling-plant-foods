% This script loads the 'clutterData_*.mat' structure created by the 'build_tables.m' script, and generates histograms showing the 
% distribution of the differences in receptor excitation values for each pixel across all images of the user-specified object category 
% for all possible dichromatic receptor comparisons. The median value for the selected category of plant food is plotted as a vertical 
% yellow line with the value printed on top of it. Also displayed is the percentage of image pixels with values falling below this median 
% value. The user needs to specify 'animal' and 'objCat' and also to uncomment the single line in lines 47-51 that corresponds to the 
% user-specified 'objCat'.

%% USER ADJUSTABLE PARAMETERS
animal = 'generalTetra'; % 'bird' or 'generalTetra'
objCat = 'immature fleshy fruit'; % 'leaf buds' or 'immature upper' or 'immature lower' or 'immature fleshy fruit' or 'immature non-fleshy fruit'
%%
close all
load(['clutterData_', animal, '.mat']); % name of file to load
clearvars -except animal objCat clutterData

if strcmp(animal, 'bird') 
    if strcmp(objCat, 'leaf buds')
        textX = [0.1 0.14 0.08 0.18 0.17 0.17]; 
    elseif strcmp(objCat, 'immature upper')
        textX = [0.15 0.1 0.1 0.07 0.19 0.14]; 
    elseif strcmp(objCat, 'immature lower')    
        textX = [0.08 0.08 0.08 0.18 0.08 0.08];
    elseif strcmp(objCat, 'immature fleshy fruit')    
        textX = [0.16 0.1 0.09 0.18 0.14 0.15];
    elseif strcmp(objCat, 'immature non-fleshy fruit')    
        textX = [0.13 0.14 0.1 0.18 0.16 0.14]; 
    end
elseif strcmp(animal, 'generalTetra')
    if strcmp(objCat, 'leaf buds')
        textX = [0.12 0.11 0.12 0.17 0.18 0.18];
    elseif strcmp(objCat, 'immature upper')
        textX = [0.17 0.11 0.08 0.18 0.17 0.15];
    elseif strcmp(objCat, 'immature lower')    
        textX = [0.09 0.08 0.08 0.18 0.09 0.09]; 
    elseif strcmp(objCat, 'immature fleshy fruit')    
        textX = [0.16 0.1 0.09 0.06 0.15 0.14]; 
    elseif strcmp(objCat, 'immature non-fleshy fruit')    
        textX = [0.14 0.13 0.09 0.18 0.18 0.17]; 
    end
end

iter = 1;
for i = 1:numel(clutterData)
        
    for j = 1:numel(clutterData(i).object)
        o = clutterData(i).object(j).objCat;
        % objects to filter analysis by (uncomment appropriate line)
        if strcmp(o, 'leaf buds') || strcmp(o, 'fiddlehead') % leaf buds
%       if strcmp(o, 'new leaves - upper') % immature upper
%       if strcmp(o, 'new leaves - lower') % immature lower
%       if strcmp(o, 'unripe fleshy fruit') || strcmp(o, 'unripe fleshy fruit capsule') % immature fleshy fruit
%       if strcmp(o, 'immature achene') || strcmp(o, 'immature nut')  || strcmp(o, 'immature seedpod') || strcmp(o, 'immature spikelet')  || strcmp(o, 'immature seed') % immature non-fleshy fruit
    
            USimgs(:,iter) = clutterData(i).contrastImgs.US;
            USfood(iter,1) = clutterData(i).object(j).USobjdiff;
            
            UMimgs(:,iter) = clutterData(i).contrastImgs.UM;
            UMfood(iter,1) = clutterData(i).object(j).UMobjdiff;
            
            ULimgs(:,iter) = clutterData(i).contrastImgs.UL;
            ULfood(iter,1) = clutterData(i).object(j).ULobjdiff;
            
            MLimgs(:,iter) = clutterData(i).contrastImgs.ML;
            MLfood(iter,1) = clutterData(i).object(j).MLobjdiff;
            
            SMimgs(:,iter) = clutterData(i).contrastImgs.SM;
            SMfood(iter,1) = clutterData(i).object(j).SMobjdiff;
            
            SLimgs(:,iter) = clutterData(i).contrastImgs.SL;
            SLfood(iter,1) = clutterData(i).object(j).SLobjdiff;
            
            iter = iter + 1;
        end
    end
end

figure
set(gcf, 'Position', [1 1 1015 50])
plotHists(USimgs, USfood, 1, textX(1))
plotHists(UMimgs, UMfood, 2, textX(2))
plotHists(ULimgs, ULfood, 3, textX(3))
plotHists(MLimgs, MLfood, 4, textX(4))
plotHists(SMimgs, SMfood, 5, textX(5))
plotHists(SLimgs, SLfood, 6, textX(6))

function [] = plotHists(imgs, foods, plotNum, textX)
    ah = subplot(1,6,plotNum);
    currentPos = get(ah,'Position');
    newPos = currentPos - [0.025*plotNum 0 0 0];
    set(ah,'Position', newPos);
    histogram(imgs)
    hold on
    xlim([-0.44 0.3])
    set(gca,'TickDir','out');
    yticks([])
    set(gca,'xticklabel',{[]})
    medianFoods = median(foods(:,1));
    x = repmat(medianFoods, 1,2);
    yl = ylim;
    y = [0 yl(2)];
    plot(x,y,'Color', [.9 .9 0], 'LineWidth', 2)
    txt1 = num2str(round(x(1), 2));
    text(x(1)-textX, 0.85*y(2), txt1)
    [h,w,z] = size(imgs);
    perc = (sum(sum(imgs < medianFoods)) / (h * w * z)) * 100;
    txt2 = [num2str(round(perc,0)), '%'];
    xl = xlim;
    text(xl(1)+0.025, 0.2*y(2), txt2)
end