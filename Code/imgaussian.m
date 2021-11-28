function I=imgaussian(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with 
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
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

% Function is written by D.Kroon University of Twente (September 2009)


I = imgaussfilt3(I,sigma);

% 
% if(~exist('siz','var')), siz=sigma*6; end
% 
% if(sigma>0)
%     % Make 1D Gaussian kernel
%     x=-ceil(siz/2):ceil(siz/2);
%     H = exp(-(x.^2/(2*sigma^2)));
%     H = H/sum(H(:));
% 
%     % Filter each dimension with the 1D Gaussian kernels\
%     if(ndims(I)==1)
%         I=imfilter(I,H, 'same' ,'replicate');
%     elseif(ndims(I)==2)
%         Hx=reshape(H,[length(H) 1]);
%         Hy=reshape(H,[1 length(H)]);
%         I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
%     elseif(ndims(I)==3)
%         if(size(I,3)<4) % Detect if 3D or color image
%             Hx=reshape(H,[length(H) 1]);
%             Hy=reshape(H,[1 length(H)]);
%             for k=1:size(I,3)
%                 I(:,:,k)=imfilter(imfilter(I(:,:,k),Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
%             end
%         else
%             Hx=reshape(H,[length(H) 1 1]);
%             Hy=reshape(H,[1 length(H) 1]);
%             Hz=reshape(H,[1 1 length(H)]);
%             I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
%         end
%     else
%         error('imgaussian:input','unsupported input dimension');
%     end
% end