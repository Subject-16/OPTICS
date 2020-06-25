clc
%close all
clear
%Simulating proagation of a gaussian beam through phase screen and SLM
% al values in SI. Initially fraunhauffer diff will be assumed
N = 600;%600by600 res
L = 5e-3;% 5 mm length of side
d     = L/N; % Source plane spacing
wv   = 532e-9;% solid state laser wavelength 
k     = 2*pi/wv;
Dz    = 0.1; % propagation distance
wst0  = 0.35e-3;% beam waist
Zr    = pi*wst0^2/(wv);% rayleigh length
wstz  = wst0*sqrt(1+(Dz/Zr)^2);
Rz    = Dz*(1+(Zr/Dz)^2);
phi   = atan(Dz/Zr);

[x1,y1]  = meshgrid((-N/2:N/2-1)*d);%imp for proper centering of coords

Psibeam = wstz/wst0*exp(-(x1.^2 + y1.^2)/wstz^2).*exp(-i*k*(x1.^2 + y1.^2)/(2*Rz))*exp(i*(phi+Dz));
figure

image(abs(Psibeam),'CDataMapping','scaled')% object plane

[Uout,x2,y2] = fresnel_prop(Psibeam,wv,d,Dz);
figure;
Inte = abs(Uout).^2;% intensity
image(Inte,'CDataMapping','scaled')% image plane

[z ] = roughness2(N,N);
Uout = Uout.*exp(1i*k*(z));
figure;
Inte = abs(Uout).^2;% intensity
image(Inte,'CDataMapping','scaled')% image plane

%again propagating the beam same distance

[Uout,x2,y2] = fresnel_prop(Uout,wv,d,Dz);
figure;
Inte = abs(Uout).^2;% intensity
image(Inte,'CDataMapping','scaled')% image plane


