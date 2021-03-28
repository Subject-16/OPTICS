% test file will expand later
clc
close all
clear

L   = 5e-3;
N   = 600;
dx  = (L/N);
ll  = 1;
g   = (-N/2:N/2-1)*dx;
Z   = 10;
lamb  = 627e-9;
alph  = 0:0.2:1.6;

[x,y] = meshgrid(g);
TT    = mod(ll*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);
del   = 0.25:0.25:4;
w     = 0.02e-3; w1 = sqrt(2)*w;
jj    = 4;
for ii = 1:length(alph)
    
Vbeam = sin(alph(ii))*exp(-(x.^2 + y.^2)/w^2) + cos(alph(ii))*exp(-(x.^2 +...
   y.^2)/w1^2).*T.*exp(1i*del(jj));
% 
% figure
% imagesc(((abs(Vbeam)).^2),'CDataMapping','scaled')
% title('Test pol');

[Uout,x2,y2] = fraunhofer_prop(Vbeam,lamb,dx,Z);
MX = max(max((abs(Uout)).^2));

subplot(3,3,ii);
imagesc(((abs(Uout)).^2)/MX,'CDataMapping','scaled')
title(alph(ii));

end

ii    = 4;
figure
for jj = 1:length(del)
    
Vbeam = sin(alph(ii))*exp(-(x.^2 + y.^2)/w^2) + cos(alph(ii))*exp(-(x.^2 +...
   y.^2)/w1^2).*T.*exp(1i*del(jj));
% 
% figure
% imagesc(((abs(Vbeam)).^2),'CDataMapping','scaled')
% title('Test pol');

[Uout,x2,y2] = fraunhofer_prop(Vbeam,lamb,dx,Z);
MX = max(max((abs(Uout)).^2));


subplot(4,4,jj);
imagesc(((abs(Uout)).^2)/MX,'CDataMapping','scaled')
title(del(jj));

end
