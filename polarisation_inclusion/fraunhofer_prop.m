function [Uout,x2,y2] = fraunhofer_prop(Uin, wvl, d3, Dz)   


if sum(isnan(Uin(:))) > 0
    disp('NaN Values in Matrix. Replacing Now');
    Uin = replace_nan1(Uin);
end
    
    

N       = size(Uin, 1);%assume square grid
k       = 2*pi/ wvl;% optical wavevector
fX      =(-N/2 : N/2-1) / (N*d3);% observation-plane coordinates
[x2,y2] = meshgrid(wvl*Dz*fX);

clear('fX');

Uout    = exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/(1i*wvl*Dz).*fftshift(fft2(Uin))*d3^2;%for 2d cases

end

