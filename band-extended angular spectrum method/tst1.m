clc
close all
clear

%% HEIST SET UP
wv    = 627e-9;% solid state laser wavelength must be orders less than d
N     = 600;
L     = 5e-3;
d     = L/(2*N);
Dz    = .05; % propagation distance
wst0  = 0.35e-3;% beam waist
E0    = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N    = 600;%600by600 res
% L    = 5e-3;% 5 mm length of side
% d    = L/N; % Source plane spacing
% wv   = 532e-9;% solid state laser wavelength 
% Dz   = 0.1; % propagation distance
% wst0 = 0.35e-3;% beam waist
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X       = (-N/2:N/2-1)*d;% defining the real space coordinates
[x,y]  = meshgrid(X,X);


Psibeam = E0*exp(-(x.^2 + y.^2)/(wst0^2));
size(Psibeam)
N       = length(Psibeam);

figure
imagesc(abs(Psibeam),'CDataMapping','scaled')
title('Input Spectrum before padding')


%  Padding the gaussian with zeros
Psibeam1 = padarray(Psibeam,N/2);
Psibeam1 = padarray(Psibeam1',N/2);
size(Psibeam1)
n       = length(Psibeam1);

X1      = (-n/2:n/2-1)*d;% defining the real space coordinates
[x1,y1] = meshgrid(X1,X1);
k       = (-n/2:n/2-1)/(n*d);
[kx,ky] = meshgrid(k);

figure
imagesc(abs(Psibeam1),'CDataMapping','scaled')
title('Input Spectrum after padding')


zc = n*d^2/(wv)*sqrt(1-(wv/(2*d))^2);
figure
plot(X1,Psibeam1(n/2,:))

close all

%% Analytical solution without padding
kv   = 2*pi/wv;
gam  = wv*Dz/(pi*wst0^2);
wz   = wst0*sqrt(1+gam^2);
Rz   = Dz*(1+1/gam^2);
phi = kv*Dz + kv/(2*Rz)*(x1.^2 + y1.^2);
Uan  = E0/(1+1i*gam)*exp(-(x1.^2+y1.^2)/(wz^2)).*exp(1i*phi);


%% HEIST DAY ANGULAR SPECTRUM METHOD(normal fft) with PADDING

[U]   = Angular_spect(Psibeam1,d,wv,Dz);

figure
imagesc(abs(U),'CDataMapping','scaled')
title('Normal Angular Spectrum z=0.1m')

figure
hold
plot(X1,abs(U(n/2,:)))
plot(X1,abs(Uan(n/2,:)))
title('Normal Angular Spectrum z=0.1m')
hold

%% Analytical solution without padding
kv   = 2*pi/wv;
gam  = wv*Dz/(pi*wst0^2);
wz   = wst0*sqrt(1+gam^2);
Rz   = Dz*(1+1/gam^2);
phi = kv*Dz + kv/(2*Rz)*(x1.^2 + y1.^2);
Uan  = E0/(1+1i*gam)*exp(-(x1.^2+y1.^2)/(wz^2)).*exp(1i*phi);

%% HEIST DAY ANGULAR SPECTRUM METHOD(NUFFT) WITH PADDING

[U]   = Band_extended_angular_spect(Psibeam1,x1,y1,d,wv,Dz);

figure
imagesc(abs(U),'CDataMapping','scaled')
title('Modified Angular Spectrum z=0.1m')

figure
hold
plot(X1,abs(U(n/2,:)))
plot(X1,abs(Uan(n/2,:)))
title('Modified Angular Spectrum z=0.1m')
hold

% %%  %%% multiplying with 2d rough phase screen
% % 
% % z  = randn(N,N);
% % U1  = U.*exp(1i*2*pi/wv*(z));
% % figure;
% % Inte = abs(U1);% intensity
% % image(Inte,'CDataMapping','scaled')% image plane
% % title('After multiplying with phase screen')
% % %%% propagating the 1m
% % Dz2 = 1;
% % 
% [U2]   = Angular_spect(U,d,wv,Dz);
% 
% figure
% imagesc(abs(U2),'CDataMapping','scaled')
% title('After propagation')
% 
% %%  generating the OAM modes
% Inte = OAM_GEN(U,d,L,10);
% 
% %% checking correlation
% [Corr,F] = Rough_2d(Inte,n);