% This program iteratively optimize the phase profile encoded on the SLM
% using Fraunhoffer diffraction under the assumption that SLM is precisely
% put on the pupil of the imaging lens.(In practice this may needs further
% approximation by adding a lens phase pattern.) The algorithm is based on
% the Gerchberg-Saxton approach. The hybrid-input-output algorithm can be
% tested as well, however, the results are not close to be convincing.

clear all;
close all;
clc;

%% Initialization of imaging parameters
% N_iter=1;
% wavefront
f = 150e-3;
lambda = 550e-9; 
k = 2*pi/lambda; % optical wavevector

% Pupil plane, use SLM coordinates
deltaS = 8e-6;
NS = 1920;
dS = (-NS/2:NS/2-1)*deltaS;
[xS yS] = meshgrid(dS);
rS = sqrt(xS.^2 + yS.^2);
% rS = rS(960-540:960+539,:);

% target normalized PSF
a = 8e8; % Gaussian function: exp(-ax^2)
r0 = 25*deltaS;
A = 1;

% for type = 0:1:5
type = 0;
switch(type)
%% Target pattern
%% Shifted Dot
% If0 = zeros(size(rS));
% If0(961,961+shiftN*10) = 1;

%% Spiral
case 0 
theta = 0:pi/100:6*pi;
R = theta.^2;

x = R.*cos(theta)/(36*pi^2);
y = R.*sin(theta)/(36*pi^2);

If0 = zeros(1920,1920);
x_img = floor(50*x);
y_img = floor(50*y);

thr = 18*deltaS;
ringMask = ones(size(rS));
ringMask(rS<thr) = 0;

for idx = 1:601
    If0(uint16(y_img(idx)+960):uint16(y_img(idx)+962),uint16(x_img(idx)+960):uint16(x_img(idx)+962)) = 1;
%     If0(uint16(y_img(idx)+961), uint16(x_img(idx)+961)) = 1;
end
If0 = If0.*ringMask;
phaseVideoName= sprintf('./results/videos/phaseVideo_Spiral.avi');

%% Symmetric dots
case 1
If0 = zeros(size(rS));
If0(959-26,959) = 1;
If0(959,959+26) = 1;
If0(959+26,959) = 1;
If0(959,959-26) = 1;
If0(959+13,959+13) = 1;
If0(959-13,959-13) = 1;
If0(959-13,959+13) = 1;
If0(959+13,959-13) = 1;
phaseVideoName= sprintf('./results/videos/phaseVideo_dots.avi');

%% Gaussian distributed ring
case 6
If0 = A/(lambda*f) * exp(-a*(rS-r0).^2); 
mask = ones(size(If0));
hWidth = 2.5*deltaS;
mask((rS-r0-hWidth) >0) = 0;
mask((rS-r0+hWidth) <0) = 0;
If0 = If0 .* mask;

%% Gaussian distributed SHIFTED ring
case 2
rS = sqrt((xS-deltaS*50).^2 + yS.^2);
If0 = A/(lambda*f) * exp(-a*(rS-r0).^2); 
mask = ones(size(If0));
hWidth = 2.5*deltaS;
mask((rS-r0-hWidth) >0) = 0;
mask((rS-r0+hWidth) <0) = 0;
If0 = If0 .* mask;
phaseVideoName= sprintf('./results/videos/phaseVideo_shiftedRing.avi');
%% Peak ring
% hWidth = 0.5*deltaS;
% If0 = ones(size(rS));
% If0((rS-r0-hWidth)> 0) = 0;
% If0((rS-r0+hWidth)< 0) = 0;
% hWidth = 2.5*deltaS;
% mask = ones(size(If0));
% mask((rS-r0-hWidth) >0) = 0;
% mask((rS-r0+hWidth) <0) = 0;
% If0 = If0 .* mask;

%% Letter E
case 3
If0 = double(imread('letterE.png'));
mask = If0;
mask = mask/max(max(mask));
phaseVideoName= sprintf('./results/videos/phaseVideo_E.avi');

%% SIGGRAPH logo
% If0 = double(imread('acmsigPSF.png'));
% mask = If0;
% mask = mask/max(max(mask));

%% Gaussian distributed double rings
case 4
x0 = 50*deltaS; y0 = 0;
If0 = A/(lambda*f) * (exp(-a*(sqrt((xS-x0).^2+(yS-y0).^2)-x0).^2) + exp(-a*(sqrt((xS+x0).^2+(yS-y0).^2)-x0).^2));
mask0 = ones(size(If0));
mask1 = ones(size(If0));
hWidth = 3*deltaS;
mask0((sqrt((xS-x0).^2 + (yS-y0).^2) -x0 - hWidth)>0) = 0;
mask0((sqrt((xS-x0).^2 + (yS-y0).^2) - x0 + hWidth)<0) = 0;
mask1((sqrt((xS+x0).^2 + (yS-y0).^2) -x0 - hWidth)>0) = 0;
mask1((sqrt((xS+x0).^2 + (yS-y0).^2) -x0 + hWidth)<0) = 0;
mask = mask1 + mask0; mask = mask/max(max(mask));

ringMask = ones(size(rS));
ringMask(rS<18*deltaS) = 0;
If0 = If0.*ringMask;
phaseVideoName= sprintf('./results/videos/phaseVideo_doubleRings.avi');
%% Gaussian distributed half rings
case 5
x0 = 50*deltaS; y0 = 0;
r0 = 50*deltaS;
yS0 = yS;
yS0(yS<=0) = 2*r0;
yS1 = yS;
yS1(yS>=0) = -2*r0;

ringMask = ones(size(rS));
ringMask(rS<18*deltaS) = 0;

If0 = A/(lambda*f) * (exp(-a*(sqrt((xS-x0).^2+(yS0-y0).^2)-r0).^2) + exp(-a*(sqrt((xS+x0).^2+(yS1-y0).^2)-r0).^2));
If0 = If0.*ringMask;
phaseVideoName= sprintf('./results/videos/phaseVideo_halfRings.avi');
end % switch type
%% Blurred PSF
h = fspecial('disk',3);
h = fspecial('gaussian',15, 1.5);
If0 = imfilter(If0, h, 'replicate');

%% Mask for hybrid-input-output approach
mask = If0; % mask for hybrid-input-output approach
mask = mask/max(max(mask));

%% Normalization of PSF
If0 = If0/sum(sum(If0)); % normalize point spread function

%% initialize video recording
% phaseUncompVideo = VideoWriter(phaseVideoName, 'Uncompressed AVI');
% psfUncompVideo = VideoWriter('./results/videos/psfVideo_gaussianRingSmth_60Hz.avi', 'Uncompressed AVI');
% phaseUncompVideo.FrameRate = 60;
% psfUncompVideo.FrameRate = 60;
% open(phaseUncompVideo);
% open(psfUncompVideo);
% for frame = 1:60

%% change iteration number
% for N_iter = 1:10
N_iter = 10;
%% Phase-retrieval optimization
% init Uf, the target optical wave distribution on image plane
Uf0 = sqrt(If0);
% fUf = rand(size(If0)) + j * rand(size(If0)); % random wave function
fUf = Uf0; 

mask_slm_hd = zeros(size(fUf)); mask_slm_hd(960-540:960+539,:) = 1; % slm frame mask
constraints = []; % constraint error
se = []; % square error
P = rand(size(If0)) .* exp(j * rand(size(If0)) * 6*pi); % random phase

for iter = 1:N_iter
%% set the ampitude of pupil function to be homogeneous
P = exp(j*angle(P));
P = P.*mask_slm_hd;

%% forward propagation
fUf = exp(j*k*rS.^2/(2*f))/(j*lambda*f) .* fftshift(fft2(fftshift(P))) * deltaS^2;
If = fUf .* conj(fUf);
If = If/sum(sum(If)); % normalization

%% image plane constraints
constraints = [constraints, sum(sum(If - If .* mask))];

%% apply image plane constraints(hybrid-input-output approach)
% beta = 0.9; 
% fUf = fUf-fUf.*(1-mask)*beta;

%% apply image plane constraints by substituting PSF modulus
fUf = exp(j*angle(fUf)) .* Uf0;

%% backward propagation
p = fUf * (j*lambda*f) ./ exp(j*k*rS.^2);
P = 6*ifftshift(ifft2(ifftshift(p))) * (deltaS * NS)^2;

%% estimate square error
se = [se,sum(sum((If-If0).^2))];
end

% additional penalty against the non-smoothness using Laplatian filtered image
% w = fspecial('laplacian', 0);
% N_pnt_iter = 100000;
% lambda_smth = 2e-4;
% mask_fUf = zeros(size(If0));
% mask_fUf(Uf0>0) = 1;
% for iter = 1:1:N_pnt_iter
% %% set the ampitude of pupil function to be homogeneous
% P = exp(j*angle(P));
% P = P.*mask_slm_hd;
% % figure;imshow(abs(P));
% %% forward propagation
% fUf = exp(j*k*rS.^2/(2*f))/(j*lambda*f) .* fftshift(fft2(fftshift(P))) * deltaS^2;
% If = fUf .* conj(fUf);
% If = If/sum(sum(If)); % normalization
% 
% %% apply image plane constraints by substituting PSF modulus
% fUg = imfilter(abs(fUf), w, 'replicate');
% fUf = exp(j*angle(fUf)) .* abs(Uf0 + lambda_smth*fUg) .* mask_fUf;
% % figure;imshow(abs(fUf));
% %% backward propagation
% p = fUf * (j*lambda*f) ./ exp(j*k*rS.^2);
% P = ifftshift(ifft2(ifftshift(p))) * (deltaS * NS)^2;
% 
% end

%% Reoptimize with image plane constraints
% N_iter_2 = 10;
% for iter = 1:N_iter_2
% %% set the ampitude of pupil function to be homogeneous
% P = P./abs(P);
% P = P.*mask_slm_hd;
% 
% %% forward propagation
% fUf = exp(j*k*rS.^2/(2*f))/(j*lambda*f) .* fftshift(fft2(fftshift(P))) * deltaS^2;
% If = fUf .* conj(fUf);
% If = If/sum(sum(If)); % normalization
% 
% %% image plane constraints
% constraints = [constraints, sum(sum(If - If .* mask))];
% 
% %% apply image plane constraints(hybrid-input-output approach)
% % beta = 0.9; 
% % fUf = fUf-fUf.*(1-mask)*beta;
% 
% %% apply image plane constraints by substituting PSF modulus
% fUf = fUf ./abs(fUf) .* Uf0;
% 
% %% backward propagation
% p = fUf * (j*lambda*f) ./ exp(j*k*rS.^2);
% P = ifftshift(ifft2(ifftshift(p))) * (deltaS * NS)^2;
% 
% %% estimate square error
% se = [se,sum(sum((If-If0).^2))];
% end

%% Simulate point spread function using the opimized pupil function P
P = P./abs(P); 

% hP = fspecial('gaussian',11, 2);
P_ph = angle(P);
% P_ph = imsharpen(P_ph);
% P_ph = imfilter(P_ph+pi, hP, 'replicate')-pi;
P = 1*exp(j*P_ph);

P = P.*mask_slm_hd;
fUf = exp(j*k*rS.^2/(2*f))/(j*lambda*f) .* fftshift(fft2(fftshift(P))) * deltaS^2;
If = fUf .* conj(fUf);

%% display results
% figure;imshow(If0,[]);title('Target PSF');
% figure;imshow(angle(P)+pi,[]);title('Phase profile on SLM');
figure;imshow(If,[]);title('Simulated Point Spread Function');
% figure;plot(se);title('Square Error');
% figure;plot(constraints);title('Constraints Error');
% [Min_se idx_se] = min(se)
% [Min_c idx_c] = min(constraints)
% If = If/sum(sum(If));
% figure;plot(-960:959,If0(960,:));title('cross section of target PSF');
% figure;plot(-960:959,If(960,:));title('cross section of simulated PSF');

%% write results
% Phase = angle(P)+pi; Phase = uint8(255*Phase/(2*pi)); 
% img = zeros(1920,1080); img = Phase(960-540:960+539,:);
% writeVideo(phaseUncompVideo, img);
% writeVideo(psfUncompVideo,If/max(max(If)));
% imwrite(img,strcat('slmPhase_gaussianRingSmth_iter_',num2str(N_pnt_iter),'.png'));
% imwrite(img,strcat('slmPhase_gaussianRing_iter_',num2str(N_iter),'.png'));
% imwrite(imadjust(If,[0 max(max(If))], [0 1]), strcat('slmPSF_Ring_iter_',num2str(N_iter),'.png'));
% end % optimization iteration

% end % frame
% close(phaseUncompVideo);
% close(psfUncompVideo);
% end

