%% Test script for testing OAM decomposition
clc
clear
close all

alph  = 0:0.2:pi/2;
del   = 0:0.25:4;

N     = 600; % N (must be even)is chosen  density for points in one direction (make even)
L     = 5e-3;
d     = L/(N);
g     = (-N/2:N/2-1)*d;
[x,y] = meshgrid(g);
w      = 0.04e-3;
MM    = 1;
w1    = w;
TT    = mod(MM*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);



lamb  = 627e-9;
Z     = 10;% proapgation distance SI

alph  = alph(1);
del   = del(1);

Vbeam  = (sin(alph)*exp(-(x.^2 + y.^2)/w^2) + cos(alph)*exp(-(x.^2 +...
  y.^2)/w1^2).*T.*exp(1i*del));

[Uout,~,~] = fraunhofer_prop(Vbeam,lamb,d,Z);
MX = max(max((abs(Uout)).^2));


figure
imagesc(((abs(Uout)).^2)/MX,'CDataMapping','scaled')
title('Input Spectrum before padding')

[MM,l_weight] = OAM_decomp(Uout,w1,x,y,d,Z,lamb,g);
