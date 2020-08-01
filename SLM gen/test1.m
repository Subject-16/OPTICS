clc
clear
close all


% Code for generating SLM holograms 29-07-2020 first edit
% size of OAM order beam increases with decrease in waist of input gaussian
N     = 600; % density for points in one direction (make even)
L     = 5e-3;
d     = L/(N);
g     = (-N/2:N/2-1)*d;
[x,y] = meshgrid(g);
E0   = 1;
w     = 0.2e-3;
lamb  = 627e-9;
kbeam = 2*pi/lamb;
kslm  = kbeam/10;
R     = 0.1; % radius of curvature of gaussian doughnut depends on 
p     = 1; 
ll    = 2; %order of OAM to be generated
Z     = 0.05; % proapgation distance
% in polar coordinates
%[r,t] = meshgrid(0:d:ll,0:pi/30:2*pi);


% %% doughnut gaussian wave
% r  = sqrt(x.^2 + y.^2);
% t  = atan(y./x);
% 
% DG = E_0*r./w.*exp(-(r/w).^2).*exp(-1i*k*r.^2/(2*R) - 1i*p*t) + 1.5;
% 
% % padding with xeros for input to angular spect
% % DG = padarray(DG,N/2);
% % DG = padarray(DG',N/2);
% 
% % n        = length(DG);
% % g        = (-n/2:n/2-1)*d;
% % [x,y]    = meshgrid(g);
% imagesc(abs(DG),'CDataMapping','scaled')
% 
% U  = Band_extended_angular_spect(DG,x,y,d,lamb,Z);
% 
% figure
% imagesc(abs(U),'CDataMapping','scaled')


%% generating gaussian beam

Psibeam = E0*exp(-(x.^2 + y.^2)/(w^2));
size(Psibeam)

figure
imagesc(abs(Psibeam),'CDataMapping','scaled')
title('Input Spectrum before padding')

%% SLM generation II

kx  = kslm/10;
phi = atan(y./x);

T   = 0.5*(1-cos(kx*x - ll*phi));

figure
imagesc(T,'CDataMapping','scaled')
title('SLM')
colormap(gray)

%% Multiplying SLM with gaussian

G = T.*Psibeam;

figure
imagesc(G,'CDataMapping','scaled')
title('Gaussian*SLM')

%fraunhaufer
[Uout,x2,y2] = fraunhofer_prop(G,lamb,d,Z);
figure
imagesc((abs(Uout)),'CDataMapping','scaled')

% % Angular Spectrum
% U  = Band_extended_angular_spect(T,x,y,d,lamb,Z);
% 
% figure
% imagesc(abs(U)/max(max(abs(Uout))),'CDataMapping','scaled')

























