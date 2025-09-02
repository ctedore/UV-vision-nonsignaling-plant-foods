% This script loads 'plantFoodsData.mat' and uses the image matrices within this structure to develop images
% for viewing on an sRGB display. The images developed by this script are designed to enable an objective  
% assessment by eye on a computer display, but are not meant for quantitative analysis. The script writes the
% following images:
% 1) grayscale images representing the receptor excitation of each of the avian photoreceptor classes; they are
%      named '1_bird_U.png', '2_bird_V.png', '3_bird_Su.png', '4_bird_Sv.png', '5_bird_M.png', '6_bird_L.png'. 
% 2) two false-color images with three of the images in (1) plugged into each of the channels of the RGB display; they
%      are named 'bird_LMS.png' and 'bird_LMU.png'.
% 3) opponent comparison images; they are named 'bird_M-L.png', bird_S-L.png', 'bird_S-M.png', 'bird_U-L.png',
%     'bird_U-M.png', and 'bird_U-S.png'.
% 4) one false-color image with S-U, M-L, and U-L opponent comparison images plugged into each of the channels of
%      the RGB display.
% 5) images showing which pixels were selected for color contrast calculations; each image shows the selected 
%      pixels for a single selection, with selected pixels at normal intensity and unselected pixels artificially darkened. 
%      All such images begin with 'obj' followed by the object's unique ID and categorization.
% Note that the values of over- and under-exposed pixels have been replaced with NaN, which appear as white pixels
% in the developed images. For a detailed explanation of the logic and mathematics used to create the images in 
% (1-4) above, see the Methods section of the associated manuscript. 

clear
%loads imageset
imageSetFilename = 'plantFoodsData.mat'; % filename for downloaded imagesets structure
% loads downloaded imagesets structure and renames it 'data'
loadedData = load(imageSetFilename); % loads data 
f = fieldnames(loadedData);
data = loadedData.(f{1:end}); % renames imageSets to 'data' for simplicity
clearvars -except data  % clears workspace of all variables except for 'data' 
numFilters = 6;
filterNames = {'bird_U'; 'bird_V'; 'bird_Su'; 'bird_Sv'; 'bird_M'; 'bird_L'};
removeGamma = 1;
checkSelect = 1;
alphaLev = 0.05; 
vFilters = 1;
vfilterNames = {'bird_DoubleCone'}; 
coeffFilename = 'birdFilterCoefficients.mat';
yDim = 1036;
xDim = 1392;
load(coeffFilename);
coeffSet = coeffSet(:,7);
numPhotorec = size(coeffSet,2);
allfilters=zeros(yDim,xDim,numFilters+numPhotorec);
missingData = [];
num_imagesets = length(data);
d = 'Images Adjusted for Viewing on a Digital Display (Not for Quantitative Analysis)';
mkdir(d)

for k = 1:num_imagesets
    disp(['imagset ', num2str(data(k).id)])
    
    for i = 1:numFilters
        allfilters(:,:,i) = data(k).images(i).image_matrix;
    end
    
    exposures = zeros(numFilters,1);
    for i = 1:numFilters
        exposures(i) = data(k).images(i).exposure;
    end
    mydata_exposCorr = cell(1,numFilters);
    exposCorr = 2000 ./ exposures(:,1);
    for i = 1:numFilters
        mydata_exposCorr{i} = allfilters(:,:,i) * exposCorr(i);
    end

    for i = 1:numPhotorec
        allfilters(:,:,numFilters+i) = coeffSet(1,i) * mydata_exposCorr{1} + ...
                           coeffSet(2,i) * mydata_exposCorr{2} + ...
                           coeffSet(3,i) * mydata_exposCorr{3} + ...
                           coeffSet(4,i) * mydata_exposCorr{4} + ...
                           coeffSet(5,i) * mydata_exposCorr{5} + ...
                           coeffSet(6,i) * mydata_exposCorr{6};
    end
  
    for i = 1:numFilters + numPhotorec
        mydata =  allfilters(:,:,i);
        mydata = mydata / mean(mydata(~isnan(mydata))); 
        mydata = mydata ./ (mydata + 1);
        allfilters(:,:,i) = mydata;
    end
    
    % NORMALIZE BY MAX ACROSS ALL IMAGES
    max_vals = zeros(numFilters+numPhotorec,1);
    for i = 1:numFilters+numPhotorec
        max_vals(i) = max(max(allfilters(:,:,i)));
    end
    max_for_norm = max(max_vals);
    for i = 1:numFilters+numPhotorec    
        allfilters(:,:,i) = allfilters(:,:,i) / max_for_norm;    
    end
    
    direct = [d, '/imagset_', num2str(data(k).id), '/'];
    mkdir(direct)
    
    % OPPONNENT COMPARISON IMAGES 
    % grayscale
    maximum(1) = max(max(abs(allfilters(:,:,1) - allfilters(:,:,3))));
    maximum(2) = max(max(abs(allfilters(:,:,1) - allfilters(:,:,5))));
    maximum(3) = max(max(abs(allfilters(:,:,1) - allfilters(:,:,6))));
    maximum(4) = max(max(abs(allfilters(:,:,5) - allfilters(:,:,6))));
    maximum(5) = max(max(abs(allfilters(:,:,3) - allfilters(:,:,5))));
    maximum(6) = max(max(abs(allfilters(:,:,3) - allfilters(:,:,6))));
    factor = 1 / max(maximum);
    [~] = opponentGray(allfilters(:,:,1), allfilters(:,:,3), factor, 1, direct, 'bird_U-S.png');
    [~] = opponentGray(allfilters(:,:,1), allfilters(:,:,5), factor, 1, direct, 'bird_U-M.png');
    [~] = opponentGray(allfilters(:,:,1), allfilters(:,:,6), factor, 1, direct, 'bird_U-L.png');
    [~] = opponentGray(allfilters(:,:,5), allfilters(:,:,6), factor, 1, direct, 'bird_M-L.png');
    [~] = opponentGray(allfilters(:,:,3), allfilters(:,:,5), factor, 1, direct, 'bird_S-M.png');
    [~] = opponentGray(allfilters(:,:,3), allfilters(:,:,6), factor, 1, direct, 'bird_S-L.png');
    
    % 3 opponnent comps in 1 image 
    oppComps(:,:,1) = allfilters(:,:,3) - allfilters(:,:,1);
    oppComps(:,:,2) = allfilters(:,:,5) - allfilters(:,:,6);
    oppComps(:,:,3) = allfilters(:,:,1) - allfilters(:,:,5);
    minm = min(min(min(oppComps)));
    if minm < 0
        oppComps = oppComps + abs(minm);
    elseif minm > 0
        oppComps = oppComps - minm;
    end
    maxm = max(max(max(oppComps)));
    oppComps = oppComps / maxm;
    for i = 1:3
        oppComps(:,:,i) = remGam(oppComps(:,:,i));
    end
    imwrite(oppComps, [direct, 'bird_S-U_M-L_U-M.png'])
    
    % REMOVE GAMMA SCALING OF sRGB DISPLAY SUCH THAT RECEPTOR EXCITATION VALUES SCALE LINEARLY WITH DISPLAY BRIGHTNESS
    if removeGamma == 1
        for i = 1:numFilters+numPhotorec
            allfilters(:,:,i) = remGam(allfilters(:,:,i));
        end
    end

    % SAVE MONOCHROME RECEPTOR EXCITATION IMAGES
    for i = 1:numFilters
        filename = [direct, num2str(i), '_', filterNames{i}, '.png'];
        imwrite(allfilters(:,:,i), filename);
    end
    for i = 1:numPhotorec
        filename = [direct, vfilterNames{i}, '.png'];
        imwrite(allfilters(:,:,numFilters+i), filename);
    end

    %  VARIOUS RGB COMBINATIONS OF MONOCHROME RECEPTOR EXCITATION IMAGES
    LMSu=zeros(yDim,xDim,3);
    LMSu(:,:,1)=allfilters(:,:,6);
    LMSu(:,:,2)=allfilters(:,:,5);
    LMSu(:,:,3)=allfilters(:,:,3);
    imwrite(LMSu, [direct, 'bird_LMS.png'])

    LMU=zeros(yDim,xDim,3);
    LMU(:,:,1)=allfilters(:,:,6);
    LMU(:,:,2)=allfilters(:,:,5);
    LMU(:,:,3)=allfilters(:,:,1);
    imwrite(LMU, [direct, 'bird_LMU.png'])
    
    % VIEW OBJECT SELECTIONS
    if checkSelect == 1
        selObs=size(data(k).selection_data, 1); 
        for i = 1:selObs
            f = cat(3, allfilters(:,:,6), allfilters(:,:,5), allfilters(:,:,1));
            alpha = zeros(yDim, xDim);
            if isempty(data(k).selection_data(i).selected_pixel_file)
            else
                obj = data(k).selection_data(i).selected_pixel_file;
                if sum(obj) == 0
                    disp(['Warning: object ', num2str(data(k).selection_data(i).id), ' contains no data']);
                    missingData = [missingData; data(k).id, data(k).selection_data(i).id];
                else
                    objData = reshape(obj, yDim, xDim);
                    alpha = objData + alpha;
                    alpha(alpha==0) = alphaLev;
                    f = f .* alpha;
                end
            end
            if isempty(data(k).selection_data(i).sunlit)
                illum = '';
            elseif data(k).selection_data(i).sunlit == 1
                illum = ' sun';
            elseif data(k).selection_data(i).sunlit == 0
                illum = ' shade';
            end
            imwrite(f, [direct, 'obj', num2str(data(k).selection_data(i).id), ' ', data(k).selection_data(i).object, '.png']);
            close all
        end        
        
    end
end

%% FUNCTIONS
% makeS opponent comparison images
function diffImgGray = opponentGray(chan1, chan2, factor, removeGamma, direct, fileName)
    diffImgGray = chan1 - chan2;
    diffImgGray = factor * diffImgGray * 0.5 + 0.5;
    if removeGamma == 1
         diffImgGray = remGam(diffImgGray);
    end
    if ~isempty(fileName)
        imwrite(diffImgGray, [direct, fileName])
    end
end

% removes gamma scaling of sRGB screen
function inputData = remGam(inputData)
    for x = 1:1036
        for y = 1:1392
            if inputData(x,y) <=  0.04045 
                inputData(x,y) = inputData(x,y) / 12.92;
            elseif inputData(x,y) >  0.04045
                 inputData(x,y) = ((inputData(x,y) + 0.055) / 1.055)^2.4;
            end
        end
     end
end