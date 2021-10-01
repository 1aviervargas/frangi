
%% Setup
addpath('../Code');
% Inside Code folder run in the Matlab command line
mex eig3volume.c

%% Soft mask generator
%Read simulated map from PDB low pass filtered with Chimera by a Gaussian 
%filter with variance 1.5  
m = ReadMRC('emd_12187_additional_1_filt.mrc');
mask = (m >= 0.02);
%Dilate the mask
se = strel('sphere',1);
mask = imdilate(double(mask),se);
%Low-pass filter the mask
mask = imgaussfilt3(double(mask),0.65);

%% Process EMDB map
% Read emd map to be process. This map can be downloaded from:
%https://www.emdataresource.org/EMD-12187
m = ReadMRC('emd_12187_additional_1.map');

% Filter Options
options.FrangiScaleRange=[2 3.2]; %Corresponding resolution range [2*2*0.529 3.2*2*0.529]
options.FrangiC = 1e-5; %we do not want to remove potential densities higher values removes noise but it may
%filter low density values.
options.FrangiScaleRatio = 0.2;
options.FrangiAlpha = 0.45;
options.FrangiBeta = 0.45;

[J W] = FrangiFilter3D(m,options);
J = J./max(J(:));
%We loss-pass filter the tube-like likelihood to avoid possible masking effects
J = imgaussfilt3(J,0.35);

%Enhanced map
WriteMRC(J.*m.*mask,0.529,'emd_12187_en.mrc');
%Scale map which is related with the local resolution map. Correspondence
%between scale map and local resolution map is LocalRes = 2*0.529*W; where
%0.529 is the sampling rate.
WriteMRC(W,0.529,'emd_12187_scale.mrc');
%EMDB map filtered by the same mask
WriteMRC(m.*mask,0.529,'emd_12187_mask.mrc');