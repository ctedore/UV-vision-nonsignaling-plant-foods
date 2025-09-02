% This script:
% 1) Reads in the transmittance spectra of the camera filters (files starting with 'Transmission....txt')
% 2) Reads in the spectral characteristics of all other camera components ('cameraComponentsSpectra.mat')
% 3) Calculates each channel's effective spectral sensitivity and saves them all in a file called 'realFilters.mat'. This file is later read in by the 
%      'lightcollected.m' script.
% 4) Calculates avian spectral sensitivities or generalized tetrachromat spectral sensitivities, depending on how the user sets the 'animal' variable
% 5) Finds a set of coefficients (all constrained to be >=0) that, when multiplied by the camera channel spectral sensitivities, and then summed, 
%      provides the best fit to the desired spectral sensitivities in the previous step
% 6) Outputs these coefficients in a .mat file called ‘birdFilterCoefficients.mat’ or ‘generalTetraCoefficients.mat’ depending on which animal is 
%      specified. This file is later read in by the ‘build_tables.m’ script.
% 7) Plots the desired and effective spectral sensitivities on top of each other and saves this as a jpeg file.

%% USER ADJUSTABLE VARIABLE
animal = 'bird'; % 'generalTetra' or 'bird'

%% Calculate effective spectral sensitivity of each camera channel
% load file with the spectral transmittance of the camera's IR blocking filter and CoastalOpt lens, as well as the spectral sensitivity of the JAI sensor
load('cameraComponentsSpectra.mat')
% load filter transmittance spectra and calculate effective sensitivity of each channel taking into account all camera components
files = dir('Transmission*.txt');
numfiles = length(files);
filters = zeros(size(camera,1),numfiles);
lambda = (300:991); 
for i = 1:numfiles
    fID = fopen(files(i).name);
    textscan(fID,'%f %f','HeaderLines',140);
    data = textscan(fID,'%f %f',920);
    fclose(fID);
    y = interp1(data{1,1}, data{1,2}, lambda).';
    y = y / 100;
    y = y .* camera.IRblockingFilter .* camera.CoastalOptLens .* camera.JAIsensor;     
    filters(:,i) = y;
end
save('realFilters.mat', 'filters')

%% 
outputCoeffFilename = strcat(animal, 'FilterCoefficients'); % name of output coefficients file
outputPlotName = strcat(animal, 'Filters'); % name of output plot file

%% Spectral sensitivities to simulate :
if strcmp(animal, 'bird')
    % U and V birds (only double cone is computational)
    lambMax = [370 409 449 449 504 568];
    oMid = [NaN NaN 427 452 526 588]; % oil droplet lambda_mid values you want to model
    oCut = [NaN NaN 413 433 507 562]; % oil droplet lambda_cut values you want to model    Clrs = [.64 0 0.8; 0.8 0 1; 0 0.48 .8; 0 0.6 1; 0 0.9 0; 1 0 0; 0.4 0.4 0.4]; % plot colors
    Clrs = [.64 0 0.8; 0.8 0 1; 0 0.48 .8; 0 0.6 1; 0 0.9 0; 1 0 0; 0.4 0.4 0.4];
    dblCone = 1; 
    lambMaxDbl = [568 568];
    oMidDbl = [468 NaN];
    oCutDbl = [452 NaN];
elseif strcmp(animal, 'generalTetra')
    % Generalized tetrachromat (all filters are computational)
    lambMax = [370 455 538 602];
    oMid = [NaN NaN NaN NaN]; % oil droplet lambda_mid values you want to model
    oCut = [NaN NaN NaN NaN]; % oil droplet lambda_cut values you want to model
    Clrs = [0.64 0 0.8; 0 0.48 0.8; 0 0.9 0; 1 0 0]; % plot colors
    dblCone = 0;
end

%% Optical media transmittances 
OMT_U = [0.0418	0.0436	0.0490	0.0597	0.0762	0.0972	0.1157	0.1367	0.1600	0.1856	0.2136	0.2431	0.2733	0.3031	0.3309	0.3565	0.3789	0.3978	0.4142	0.4279	0.4395	0.4491	0.4573	0.4644	0.4708	0.4774	0.4834	0.4894	0.4952	0.5007	0.5066	0.5126	0.5185	0.5247	0.5313	0.5383	0.5455	0.5525	0.5596	0.5672	0.5750	0.5831	0.5911	0.5989	0.6064	0.6139	0.6210	0.6276	0.6349	0.6417	0.6478	0.6542	0.6605	0.6666	0.6723	0.6769	0.6803	0.6836	0.6875	0.6918	0.6971	0.7011	0.7042	0.7064	0.7090	0.7120	0.7162	0.7211	0.7259	0.7298	0.7326	0.7344	0.7376	0.7411	0.7446	0.7477	0.7505	0.7531	0.7555	0.7581	0.7608	0.7634	0.7659	0.7682	0.7704	0.7725	0.7745	0.7764	0.7787	0.7809	0.7827	0.7841	0.7858	0.7875	0.7892	0.7905	0.7921	0.7935	0.7950	0.7959	0.7969	0.7979	0.7991	0.7999	0.8004	0.8013	0.8024	0.8033	0.8041	0.8047	0.8054	0.8062	0.8073	0.8083	0.8095	0.8107	0.8117	0.8129	0.8139	0.8151	0.8165	0.8180	0.8196	0.8208	0.8219	0.8232	0.8246	0.8262	0.8279	0.8299	0.8320	0.8340	0.8358	0.8376	0.8397	0.8420	0.8441	0.8461	0.8475	0.8487	0.8496	0.8504	0.8512	0.8523	0.8534	0.8544	0.8549	0.8555	0.8561	0.8567	0.8574	0.8585	0.8598	0.8606	0.8612	0.8619	0.8626	0.8635	0.8640	0.8642	0.8651	0.8668	0.8683	0.8690	0.8698	0.8707	0.8719	0.8733	0.8747	0.8760	0.8772	0.8780	0.8781	0.8784	0.8798	0.8816	0.8836	0.8848	0.8855	0.8864	0.8877	0.8888	0.8899	0.8909	0.8914	0.8918	0.8928	0.8936	0.8941	0.8942	0.8949	0.8956	0.8968	0.8983	0.8996	0.9011	0.9022	0.9030	0.9034	0.9040	0.9055	0.9061	0.9063	0.9068	0.9073	0.9078	0.9084	0.9091	0.9098	0.9103	0.9111	0.9114	0.9124	0.9136	0.9144	0.9150	0.9156	0.9167	0.9171	0.9172	0.9177	0.9184	0.9190	0.9191	0.9199	0.9208	0.9211	0.9216	0.9217	0.9228	0.9231	0.9227	0.9231	0.9238	0.9245	0.9247	0.9249	0.9256	0.9257	0.9257	0.9252	0.9260	0.9275	0.9278	0.9276	0.9279	0.9280	0.9277	0.9278	0.9283	0.9292	0.9300	0.9304	0.9303	0.9302	0.9312	0.9324	0.9331	0.9342	0.9354	0.9358	0.9356	0.9357	0.9356	0.9353	0.9354	0.9353	0.9344	0.9343	0.9341	0.9337	0.9346	0.9360	0.9365	0.9372	0.9384	0.9394	0.9396	0.9406	0.9415	0.9421	0.9428	0.9433	0.9437	0.9448	0.9455	0.9458	0.9461	0.9465	0.9459	0.9459	0.9466	0.9472	0.9468	0.9459	0.9457	0.9461	0.9462	0.9464	0.9471	0.9491	0.9500	0.9502	0.9500	0.9504	0.9512	0.9511	0.9506	0.9512	0.9515	0.9516	0.9512	0.9507	0.9504	0.9513	0.9523	0.9532	0.9537	0.9543	0.9553	0.9562	0.9567	0.9567	0.9570	0.9581	0.9586	0.9588	0.9591	0.9602	0.9613	0.9616	0.9610	0.9606	0.9608	0.9613	0.9610	0.9610	0.9614	0.9620	0.9614	0.9605	0.9602	0.9613	0.9624	0.9631	0.9634	0.9642	0.9648	0.9650	0.9648	0.9655	0.9656	0.9659	0.9661	0.9657	0.9659	0.9665	0.9667	0.9666	0.9669	0.9666	0.9662	0.9666	0.9666	0.9665	0.9665	0.9665	0.9660	0.9655	0.9658	0.9665	0.9677	0.9683	0.9690	0.9688	0.9688	0.9694	0.9694	0.9697	0.9708	0.9715	0.9708	0.9702	0.9709	0.9713	0.9721	0.9729	0.9734	0.9742	0.9746	0.9748	0.9747	0.9751	0.9756	0.9755	0.9760	0.9763	0.9766	0.9777	0.9773	0.9787].';
OMT_V = [0.0404	0.0404	0.0411	0.0428	0.0457	0.0500	0.0554	0.0624	0.0706	0.0798	0.0899	0.1006	0.1116	0.1225	0.1332	0.1433	0.1530	0.1622	0.1708	0.1790	0.1868	0.1941	0.2012	0.2082	0.2153	0.2226	0.2302	0.2379	0.2460	0.2546	0.2637	0.2732	0.2832	0.2937	0.3047	0.3165	0.3286	0.3413	0.3541	0.3675	0.3812	0.3954	0.4095	0.4236	0.4375	0.4516	0.4657	0.4801	0.4944	0.5089	0.5226	0.5351	0.5466	0.5571	0.5667	0.5762	0.5855	0.5939	0.6013	0.6078	0.6136	0.6195	0.6255	0.6313	0.6370	0.6426	0.6470	0.6510	0.6551	0.6597	0.6652	0.6709	0.6760	0.6806	0.6852	0.6895	0.6937	0.6979	0.7018	0.7060	0.7100	0.7137	0.7172	0.7206	0.7240	0.7272	0.7304	0.7334	0.7362	0.7389	0.7413	0.7435	0.7456	0.7477	0.7498	0.7516	0.7536	0.7553	0.7569	0.7584	0.7599	0.7613	0.7626	0.7639	0.7651	0.7661	0.7669	0.7676	0.7683	0.7691	0.7699	0.7708	0.7716	0.7724	0.7733	0.7743	0.7756	0.7771	0.7788	0.7807	0.7827	0.7849	0.7874	0.7902	0.7932	0.7961	0.7990	0.8019	0.8047	0.8077	0.8106	0.8133	0.8161	0.8185	0.8208	0.8231	0.8253	0.8275	0.8296	0.8318	0.8340	0.8360	0.8379	0.8396	0.8414	0.8431	0.8448	0.8465	0.8481	0.8496	0.8510	0.8522	0.8534	0.8546	0.8560	0.8573	0.8587	0.8599	0.8612	0.8627	0.8641	0.8654	0.8664	0.8675	0.8688	0.8701	0.8711	0.8718	0.8725	0.8734	0.8740	0.8747	0.8756	0.8768	0.8775	0.8778	0.8780	0.8782	0.8783	0.8788	0.8794	0.8803	0.8811	0.8824	0.8840	0.8857	0.8874	0.8890	0.8907	0.8926	0.8943	0.8957	0.8971	0.8986	0.8998	0.9006	0.9013	0.9020	0.9027	0.9034	0.9039	0.9043	0.9049	0.9055	0.9059	0.9065	0.9070	0.9076	0.9083	0.9087	0.9091	0.9096	0.9102	0.9106	0.9108	0.9112	0.9115	0.9119	0.9119	0.9119	0.9119	0.9121	0.9124	0.9126	0.9130	0.9137	0.9142	0.9148	0.9157	0.9168	0.9181	0.9193	0.9202	0.9208	0.9215	0.9222	0.9228	0.9236	0.9243	0.9249	0.9255	0.9262	0.9271	0.9283	0.9295	0.9305	0.9313	0.9322	0.9330	0.9337	0.9343	0.9348	0.9353	0.9357	0.9358	0.9357	0.9358	0.9359	0.9358	0.9358	0.9358	0.9360	0.9362	0.9362	0.9363	0.9364	0.9367	0.9370	0.9373	0.9378	0.9385	0.9392	0.9401	0.9410	0.9420	0.9430	0.9440	0.9449	0.9457	0.9462	0.9465	0.9468	0.9471	0.9472	0.9473	0.9475	0.9477	0.9478	0.9479	0.9481	0.9485	0.9491	0.9498	0.9506	0.9513	0.9520	0.9526	0.9531	0.9535	0.9541	0.9546	0.9551	0.9558	0.9564	0.9574	0.9584	0.9592	0.9600	0.9609	0.9618	0.9626	0.9635	0.9640	0.9645	0.9649	0.9650	0.9650	0.9651	0.9653	0.9654	0.9655	0.9654	0.9653	0.9653	0.9654	0.9657	0.9662	0.9669	0.9676	0.9683	0.9690	0.9695	0.9699	0.9706	0.9713	0.9718	0.9721	0.9721	0.9721	0.9724	0.9725	0.9725	0.9728	0.9733	0.9734	0.9736	0.9737	0.9738	0.9739	0.9740	0.9741	0.9745	0.9749	0.9752	0.9754	0.9758	0.9763	0.9767	0.9770	0.9774	0.9780	0.9787	0.9792	0.9798	0.9805	0.9812	0.9814	0.9813	0.9813	0.9816	0.9819	0.9819	0.9818	0.9818	0.9820	0.9821	0.9822	0.9824	0.9828	0.9834	0.9839	0.9848	0.9858	0.9867	0.9873	0.9876	0.9874	0.9874	0.9873	0.9871	0.9868	0.9865	0.9864	0.9865	0.9867	0.9874	0.9882	0.9895	0.9905	0.9911	0.9915].';
OMT_mean = (OMT_U + OMT_V) / 2;
OMT_U = vertcat(OMT_U, ones(481-length(OMT_U), 1));
OMT_mean = vertcat(OMT_mean, ones(481-length(OMT_mean), 1));
    
%%
filters = filters(1:481,:); % cuts off spectra above 780 nm where IR blocking filter effectively blocks all light
numFilters = size(filters,2);

%% Compute desired spectral sensitivity curves
figure
w = (300:780).';
vpig = zeros(length(w), length(lambMax)+dblCone);

% Single cones
for i = 1:length(lambMax)
    % opsin
    S = lambMax(i);
    a = 0.8795 + 0.0459 * exp(-((S - 300)^2) / 11940);
    b = -40.5 + 0.195 * S;
    alpha = 1./ (exp(69.7 .* (a - S./w)) + exp(28 .* (0.922 - (S ./ w))) + exp(-14.9 .* (1.104 - (S ./ w))) + 0.674);
    beta = 0.26 .* exp(-(((w - (189 + 0.315 .* S)) ./ b).^2));
    opsin = alpha + beta;
    
    % oil drop
    if isnan(oMid(i))
        vpig(:,i) = opsin;
    else
        oil = exp(-2.93 .* exp(-2.89 .* (0.5./(oMid(i)-oCut(i))) .* (w-oCut(i))));
        vpig(:,i) = opsin .* oil;
    end
    
    % OMT 
    vpig(:,i) = vpig(:,i) .* OMT_U;
end

% Double cones
if dblCone == 1
    vpigDbl = zeros(length(w),2);
    for i = 1:2
        % opsin
        S = lambMaxDbl(i);
        a = 0.8795 + 0.0459 * exp(-((S - 300)^2) / 11940);
        b = -40.5 + 0.195 * S;
        alpha = 1./ (exp(69.7 .* (a - S./w)) + exp(28 .* (0.922 - (S ./ w))) + exp(-14.9 .* (1.104 - (S ./ w))) + 0.674);
        beta = 0.26 .* exp(-(((w - (189 + 0.315 .* S)) ./ b).^2));
        opsin = alpha + beta;
        
        % oil drop
        if isnan(oMidDbl(i))
            vpigDbl(:,i) = opsin;
        else
            oil = exp(-2.93 .* exp(-2.89 .* (0.5./(oMidDbl(i)-oCutDbl(i))) .* (w-oCutDbl(i))));
            vpigDbl(:,i) = opsin .* oil;
        end
        
        % OMT 
        vpigDbl(:,i) = vpigDbl(:,i) .* OMT_mean;
    end
    vpig(:,length(lambMax)+1) = sum(vpigDbl,2);
end 

%% Compute Computational Filters
coeffSet = zeros(numfiles, length(lambMax+dblCone));
for i = 1:size(vpig,2)
    if (strcmp(animal, 'bird') && i == 7) || strcmp(animal, 'generalTetra')
        coeff = lsqlin(filters, vpig(:,i), [], [], [], [], zeros(1,numFilters), [], []);
        coeffSet(:,i) = coeff;
        synthOp = sum(filters .* coeff.', 2);
        plot(w, synthOp / sum(synthOp), ':', 'Color', Clrs(i,:), 'LineWidth', 1.5)
        hold on
    end
    plot(w, vpig(:,i) / sum(vpig(:,i)), 'Color', Clrs(i,:))
    hold on
    yticks([])
    xticks([300 400 500 600 700])
    xlim([300 750])
end

%% Combine Real and Computational Filters (birds only)
if strcmp(animal,'bird')
    realCoeffs = [ ...
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
    coeffSet(:,1:6) = realCoeffs;
    for i = 1:6
        plot(w, filters(:,i) / sum(filters(:,i)), ':', 'Color', Clrs(i,:), 'LineWidth', 1.5);
        hold on
    end
end

%% Save coefficients, plot and save filters 
set(gcf, 'Position', [1 1 525 280]) 
print([outputPlotName, '.jpg'], '-r600', '-djpeg')
save([outputCoeffFilename, '.mat'], 'coeffSet')
