close all
clear
clc

%generating random numbers
a=1;
b=0;
f=a.*randn(501,1)+b;
%finding autocorrelation funtion using Scaled FFT 
%multiplying elements by elements of fft and ifft
f2=fft(f).*ifft(f);
f1=f2.';
M=500;
m=-M/2:1:M/2;
n=-M:1:M;
m2=m.*m;
%scaling factor
d=-.02;
%creating the g function
x=1i*pi*d/M;
h=exp(-x*m2);
g=h.*f1;
b=zeros(1,M/2);
%padding symmetrically with zeros
gm=[b g b];
%creating h function
n2=n.*n;
hm=exp(x*n2);
% taking fft of g and h
G=ifft(gm);
H=ifft(hm);
%multiplying element by element and doing fft
F=fft(G.*H);
%rearrangement
F1=fftshift(F);
%multiplying by factor
h1=exp(x*m2);
h2=[b h1 b];
F2=h2.*F1;
%calculating autocorrelation function with normal FFT
yy2 = ifft(ifft(f).*fft(f));
yy3=fftshift(yy2);
%truncating F2 with M elements
F3=F2(M/2:3*(M/2));
k=d*m;

subplot(2,2,3);
plot(k,F3)
title('Finding correlation length using scaled version of Wiener Khinchin theorem')
xlabel('x in pixels')
ylabel('Autocorrelation')
hold on
%plotting straight line at first minima
x=-1.24;
y=-0.2:1;
plot(x*ones(size(y)),y)
x1=1.24;
y1=-.2:1;
plot(x1*ones(size(y1)),y1)
subplot(2,2,2);
plot(m,yy3)
title('Autocorrelation using Wiener Khinchin theorem')
xlabel('x in pixels')
ylabel('Autocorrelation')

subplot(2,2,1);
plot(m,f)
title('Gaussian random sequence')
xlabel('x in pixels')
ylabel('height')
%interpolating the surface
xq=-M/2:.05:M/2;
pp=spline(m,f,xq);

subplot(2,2,4);
plot(xq,pp)
title('interpolated surface')
xlabel('x in pixels')
ylabel('height')

%% finding the minima in the autocorrelation function to find the correlation length

[Minim,MinIdx] = findpeaks(real(-F3));
figure
plot(k,real(F3))
hold
plot(k(MinIdx),real(F3(MinIdx)),'*')
hold

Minima = k(MinIdx);

% finding the two minimas around zero Corr is the correlation length in
% pixels
for ii = 1:length(MinIdx)-1
     if Minima(ii+1)*Minima(ii)<0
         Corr = k(MinIdx(ii+1))-k(MinIdx(ii));
     end
end

abs(Corr)

