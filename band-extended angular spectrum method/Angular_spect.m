% Angular Spectrum Method of diffracction based on
%Carbajal-Domínguez, Adrian et al. “Numerical calculation of near 
%field scalar diffraction using angular spectrum of plane 
% waves theory and FFT.” (2010).

function [Phi] = Angular_spect(x,Dx,lamb,Z)
% x is the input 2d square matrix 
% Dx is the spacing in real domain
% lamb is the wavelength(SI)
% Z is the propagation distance

N  = length(x);
fc = 1/(lamb)*1/sqrt(1 + (Z/(N*Dx))^2); % theoretical fourier limit on local fourier frequency

% figure
% imagesc(x,'CDataMapping','scaled')
% title('input image')

%peforming the 2d fourier transform
A1 = fft2(fftshift(x));
A  = fftshift(A1);
k  = (-N/2:N/2-1)/(N*Dx);

% figure
% imagesc(real(A),'CDataMapping','scaled')
% title('Fourier Spectrum')

%a 1d slice of the spectrum
% figure
% plot(k,A(N/2,:))


%makinig the G,H function 

[kx,ky] = meshgrid(k,k);
p       = sqrt(1/lamb^2- kx.^2 - ky.^2);
G       = exp(1i*2*pi*Z.*p);
for ii = 1:length(kx)
    for jj = 1:length(ky)
        if(abs(kx(ii,jj))>fc || abs(ky(ii,jj))>fc)
            G(ii,jj) = 0;
        end
    end
end


H       = ifftshift(A1.*G);
Phi     = ifftshift(ifft2(H));
Phi     = (Phi);
%Phi     = Phi/max(abs(Phi));


end

