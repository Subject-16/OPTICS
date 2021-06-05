%% Test script for testing OAM decomposition
clc
clear
%close all

alph  = 0:0.2:pi/2;
del   = 0:0.25:4;

N     = 1400; % N (must be even)is chosen  density for points in one direction (make even)
L     = 5e-3;
d     = L/(N);
g     = (-N/2:N/2-1)*d;
[x,y] = meshgrid(g);
w      = 0.03e-3;
MM    = 0.5;
w1    = w;
w2    = 0.3e-3;
TT    = mod(MM*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);

l_maxima_alpha = zeros(size(alph));
l_maxima_del   = zeros(size(del));

TT    = mod(MM*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);

lamb  = 627e-9;
Z     = 10;% proapgation distance SI
for ii = 1:length(alph)

alph1  = 0.9;
del1   = del(1);

%OV + gaussian beam 
% Vbeam  = (sin(alph1)*exp(-(x.^2 + y.^2)/w2^2) + cos(alph1)*exp(-(x.^2 +...
%   y.^2)/w2^2).*T.*exp(1i*del1));

p     = 0;
r     = sqrt(x.^2 + y.^2);


% LG beam for checking decomposition
Vbeam = sqrt(2*factorial(p)./(pi*w2*gamma( abs(MM + p) + 1 ))).*(r.*sqrt(2)./w2).^(abs(MM))...
     .*exp(-((r.^2)./w2^2)).*T;
% 
% MM2   = 1; 
% TT    = mod(MM2*(atan2(y,x) + pi),2*pi);
% T     = exp(1i*TT); 
%  
% Vbeam2 = sqrt(2*factorial(p)./(pi*w2*gamma( abs(MM2 + p) + 1 ))).*(r.*sqrt(2)./w2).^(abs(MM(ii)))...
%      .*exp(-((r.^2)./w2^2)).*T; 
% 
% Vbeam = Vbeam1 + Vbeam1;
 
%Vbeam =  (Vbeam)/max(abs((Vbeam(:))));

% %%%%%%%comment this when unsing LG beam %%%%%%%%%%%%%
% [Vbeam,~,~] = fraunhofer_prop(Vbeam,lamb,d,Z);
% MX = max(max((abs(Vbeam)).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imagesc(((abs(Vbeam)).^2),'CDataMapping','scaled')
title('Input Spectrum before padding')

% can use either input
[MM,l_maxima_alpha(ii)] = OAM_decomp_LG_basis(alph1,Vbeam,w2,x,y,d,Z,lamb,g);


end
    

% for ii = 1:length(del)
% 
% alph1  = alph(4);
% del1   = del(ii);
% 
% Vbeam  = (sin(alph1)*exp(-(x.^2 + y.^2)/w^2) + cos(alph1)*exp(-(x.^2 +...
%   y.^2)/w1^2).*T.*exp(1i*del1));
% 
% [Uout,~,~] = fraunhofer_prop(Vbeam,lamb,d,Z);
% MX = max(max((abs(Uout)).^2));
% 
% [MM,l_maxima_del(ii)] = OAM_decomp_LG_basis(Uout,w2,x,y,d,Z,lamb,g);
% 
% end

% figure
% plot(alph,l_maxima_alpha,'r-*')
% figure
% plot(del,l_maxima_del,'b-*')
% 
% % 
% figure
% imagesc(((abs(Uout)).^2)/MX,'CDataMapping','scaled')
% title('Input Spectrum before padding')



