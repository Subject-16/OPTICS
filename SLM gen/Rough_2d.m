% making a 2d version of roughness profile based on Nascov et al 2009
% on running the code the varaible Corr gives the value of correlation
% length for the rough surface


function [Corr,F] = Rough_2d(f)% N = number of points , f is the input plane

%f       = randn(N,N);
%plotting the gaussian random surface
imagesc(log(abs(f)),'CDataMapping','scaled')
title('2d roughness profile ')

f2      = (fft2(f).*ifft2(f));% Wiener Kinchin Theorem for estimating autocovariance
f3      = ifft2(f2);
f2shift = fftshift(f3);
alpha   = 0.02;
M       = length(f);

% figure
% imagesc([249.9,250.1],[249.9,250.1],real(f2shift),'CDataMapping','scaled')
% title('Correlation using Weiner Khinchin Theorem(Small dot at the centre will be scaled in next figures)')

%% preparing a scaled inverse fourier transform which will be used to scale the autocovariance function
% F(n) = summed over m(f(m)*exp(-i*2*pi/M*alpha*m*n));

% creating a matrix coordinate base
n        = -M:M-1;
m        = -M/2:M/2-1;
kscaled  = alpha*m;
[xn,yn]  = meshgrid(n,n);
[xm,ym]  = meshgrid(m,m);

w  = exp(1i*2*pi/M);
h  = exp(1i*2*pi/M.*(alpha*(xn.^2 + yn.^2)/2));
g  = exp(-1i*2*pi/M.*(alpha/2*(xm.^2 + ym.^2))).*f2;

% padding zeros to expand the g matrix
b1       = zeros(M/2,M);
b2       = zeros(M/2,M/2);
g        = [b2,b1,b2;b1',g,b1';b2,b1,b2];
f4       = fftshift(fft2(ifft2(g).*ifft2(h)));

f5       = exp(-1i*2*pi/M.*((xm.^2 + ym.^2)*alpha/2));
f5       = [b2,b1,b2;b1',f5,b1';b2,b1,b2];

F1       = f4.*f5;
F        = F1(M/2:3*(M/2)-1,M/2:3*(M/2)-1);
F2       = real(F(M/2,:));
%kscaled  = alpha*m;

figure
imagesc(real(F),'CDataMapping','scaled')% object plan
colorbar
title('Correlation using Scaled Fourier Transform alpha=0.02 ')



% plotting a 1 profile for better view
figure
plot(kscaled,F2)
xlabel('pixels')
ylabel('Correlation')
title('1 profile of the Scaled Correlation')
%% finding the minima in the autocorrelation function to find the correlation length


[Minim,MinIdx] = findpeaks(-F2);
hold
plot(kscaled(MinIdx),F2(MinIdx),'*')
hold

Minima = kscaled(MinIdx);

% finding the two minimas around zero Corr is the correlation length in
% pixels
for ii = 1:length(MinIdx)-1
     if Minima(ii+1)*Minima(ii)<0
         Corr = kscaled(MinIdx(ii+1))-kscaled(MinIdx(ii))
     end
end


end

