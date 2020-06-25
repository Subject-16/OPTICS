function [Uout,x2,y2] = fresnel_prop(Uin, wvl, d3, Dz)   


N       = size(Uin, 1);%assume square grid
k       = 2*pi/ wvl;% optical wavevector
fX      =(-N/2 : N/2-1) ;% observation-plane coordinates
[x2,y2] = meshgrid(wvl*Dz*fX/(N*d3));
[x1,y1] = meshgrid(fX*d3);
clear('fX');

Uout    = exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/(1i*wvl*Dz).*fftshift(fft2(Uin.*exp(1i*k/(2*Dz)*(x1.^2 + y1.^2))))*d3^2; %for 2d cases

end