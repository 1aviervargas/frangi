function [Iout,Scale]=FrangiFilter3D_ori(I,options,Mask)
%
% This function FRANGIFILTER3D uses the eigenvectors of the Hessian to
% compute the likeliness of a cryo-EM map region to tube-like structures.
%
% [J,Scale,Vx,Vy,Vz] = FrangiFilter3D(I, Options)
%
% Inputs,
%   I : The input cryo-EM map
%   Options : Struct with input options,
%       .FrangiScaleRange : The range of sigmas used, default [1 2] in px
%                           units
%       .FrangiScaleRatio : Step size between sigmas, default 0.1
%       .FrangiAlpha : Frangi vesselness constant, treshold on Lambda2/Lambda3
%					   determines if its a line(vessel) or a plane like structure
%					   default .45;
%       .FrangiBeta  : Frangi vesselness constant, which determines the deviation
%					   from a blob like structure, default .45;
%       .FrangiC     : If a solvent Mask is not provided (third input optional parameter): 
%					   Frangi vesselness constant that gives the threshold between 
%					   eigenvalues of noise and vessel structure. A thumb rule is dividing the 
%					   the greyvalues of the vessels by 4 till 6, default
%					   0.1. If a solvent masl is provided: quantile of the
%					   empiral distrubution of R_{kappa} calculated using
%					   volxels outsie de mask only. A Good value is
%					   0.9
%       .verbose : Show debug information, default true
%   Mask: Optional parameter which corresponds to a binary solvent mask
%   where 1 means protein and 0 background.
%
% Outputs,
%   J : The tubular enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Vx,Vy,Vz: Matrices with the direction of the smallest eigenvector, pointing
%				in the direction of the line/vessel.
%
% Literature, 
%	Manniesing et al. "Multiscale Vessel Enhancing Diffusion in 
%		CT Angiography Noise Filtering"
%
% Example,
%   % compile needed mex file
%   mex eig3volume.c
%
%   load('ExampleVolumeStent');
%   
%   % Frangi Filter the stent volume
%   options.FrangiScaleRange=[1 1];
%   Vfiltered=FrangiFilter3D(V,options);
%
%   % Show maximum intensity plots of input and result
%   figure, 
%   subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
%   subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
%   subplot(2,2,3), imshow(V(:,:,100),[])
%   subplot(2,2,4), imshow(Vfiltered(:,:,100),[])
%
%   Javier Vargas, 01/10/21
%   jvargas@ucm.es
%   Copyright 2021, Universidad Complutense de Madrid
%   $ Revision: 1.0.0.0
%   $ Date: 01/10/21
%
%Copyright 2021 Javier Vargas @UCM
%
%Redistribution and use in source and binary forms, with or without modification, 
%are permitted provided that the following conditions are met:
%
%1. Redistributions of source code must retain the above copyright notice, 
%this list of conditions and the following disclaimer.
%
%2. Redistributions in binary form must reproduce the above copyright notice, 
%this list of conditions and the following disclaimer in the documentation 
%and/or other materials provided with the distribution.
%
%3. Neither the name of the copyright holder nor the names of its contributors 
%may be used to endorse or promote products derived from this software without 
%specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
%USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Written by D.Kroon University of Twente (May 2009)

defaultoptions = struct('FrangiScaleRange', [1 2], 'FrangiScaleRatio', 0.1, 'FrangiAlpha', 0.45, 'FrangiBeta', 0.45, 'FrangiC', 0.1, 'verbose',true);

% Process inputs
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('FrangiFilter3D:unknownoption','unknown options found');
    end
end

% Use single or double for calculations
if(~isa(I,'double')), I=single(I); end

sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);
sigmas = sort(sigmas, 'ascend');

% Frangi filter for all sigmas
for i = 1:length(sigmas),
    % Show progress
    if(options.verbose)
        disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigmas(i));

    if(sigmas(i)>0)
        % Correct for scaling
        c=(sigmas(i)^0.5);
        Dxx = c*Dxx; Dxy = c*Dxy;
        Dxz = c*Dxz; Dyy = c*Dyy;
        Dyz = c*Dyz; Dzz = c*Dzz;
    end
    
    [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    
    % Free memory
    clear Dxx Dyy  Dzz Dxy  Dxz Dyz;

    % Calculate absolute values of eigen values
    LambdaAbs1=abs(Lambda1);
    LambdaAbs2=abs(Lambda2);
    LambdaAbs3=abs(Lambda3);

    % The Vesselness Features
    Ra=LambdaAbs2./LambdaAbs3;
    Rb=LambdaAbs1./sqrt(LambdaAbs2.*LambdaAbs3);

    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D
    S = sqrt(LambdaAbs1.^2+LambdaAbs2.^2+LambdaAbs3.^2);
    A = 2*options.FrangiAlpha^2; B = 2*options.FrangiBeta^2;
	    
    if(nargin>2)
        q = quantile(S(~(Mask>0)),options.FrangiC);
        C = 2*(q)^2;
    else        
        C = 2*options.FrangiC^2;
    end
    
    % Free memory
    clear LambdaAbs1 LambdaAbs2 LambdaAbs3

    %Compute Vesselness function
    expRa = (1-exp(-(Ra.^2./A)));
    expRb =    exp(-(Rb.^2./B));
    expS  = (1-exp(-S.^2./C));
    %keyboard
    % Free memory
    clear S A B C Ra Rb

    %Compute Vesselness function
    Voxel_data = expRa.* expRb.* expS;
    
    % Free memory
    clear expRa expRb expRc;
    
    Voxel_data(Lambda2 > 0)=0; Voxel_data(Lambda3 > 0)=0;
        
    % Remove NaN values
    Voxel_data(~isfinite(Voxel_data))=0;

    % Add result of this scale to output
    if(i==1)
        Iout=Voxel_data;
        if(nargout>1)
            Scale = 2*max(sigmas)*ones(size(I),class(Iout));
        end

    else
        if(nargout>1)
            Scale(Voxel_data>Iout)=2*sigmas(i);            
        end        % Keep maximum filter response

        Iout=max(Iout,Voxel_data);
    end
end
