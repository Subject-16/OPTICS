% code for comparing analytical and numerical estimate to fraunhoffer
% diffration ALL SI units
close all
clear all

N     = 512; % no. of grid points
L     = 7.5e-3; % total length of side is SI
d     = L/N; % Source plane spacing
D     = 1e-3; % diameter of grid spacing
wv1   = 1e-6;% wavelength 
k     = 2*pi/wv1;
Dz    = 10; % propagation distance 

[x1,y1]  = meshgrid((-N/2:N/2-1)*d);%imp for proper centering of coords
% SINGLE SLIT
% Urect = rectslit(x1,y1,12*d);
% figure;
% image(Urect,'CDataMapping','scaled')% source plane
% [Uout,x2,y2] = fraunhofer_prop(Urect,wv1,d,Dz);
% figure;
% Inte = abs(Uout).^2;% intensity
% %image(Inte,'CDataMapping','scaled')% image plane
% %plot(x2(y2==0),Inte(y2==0))
%CIRCULAR SLIT
Ucirc    = circ(x1,y1,D);%input field
figure;
image(Ucirc,'CDataMapping','scaled')% source plane
[Uout,x2,y2] = fraunhofer_prop(Ucirc,wv1,d,Dz);
figure;
Inte = abs(Uout).^2;% intensity
image(Inte,'CDataMapping','scaled')% image plane
yslc = 0;%slice of y needed
Inte = abs(Uout(y2==yslc)).^2;
plot(x2(1,:),Inte)
% 



%analytical result

%Uout_th = exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/ (1i*wv1*Dz)*D^2*pi/4.*jinc(D*sqrt(x2.^2+y2.^2)/(wv1*Dz));


