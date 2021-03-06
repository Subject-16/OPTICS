% clc
% clear
% close all

% SLM generation Based on Thesis on 'Applications of Laguerre-Gaussian beams
% and Bessel beams to both nonlinear optics and atom optics' 
% by Jochen Arlt

function [x,Uout,xx,UU] = test1(n,ll,uu,alpha,phi)
%% Code for generating SLM holograms 29-07-2020 first edit
% size of OAM order beam increases with decrease in waist of input gaussian

% if Flag is 0 less rough screen will be used if flag is 1, roughness scale
% will be higher
flag  = 0;
%n     = 10; %must be even || roughness number
Nn    = n^2;
N     = 600; % N (must be even)is chosen  density for points in one direction (make even)
L     = 5e-3;
d     = L/(N);
g     = (-N/2:N/2-1)*d;
[x,y] = meshgrid(g);
E0    = 100;
w     = .35e-3;
lamb  = 627e-9;
kbeam = 2*pi/lamb;
kslm  = kbeam/10;
R     = 0.1; % radius of curvature of gaussian doughnut depends on 
p     = 1; 
% ll    = 1; %order of OAM to be generated
Z     = 10;% proapgation distance SI

% in polar coordinates
%[r,t] = meshgrid(0:d:ll,0:pi/30:2*pi);



%% generating gaussian beam

Psibeam = E0*exp(-(x.^2 + y.^2)/(w^2));
% 
% figure
% imagesc(abs(Psibeam),'CDataMapping','scaled')
% title('Input Spectrum before padding')

%% Multiplying with Gaussian Roughness screen
% inherently randn produces random nos with 0 mean and 1 variance
% Making a less rough matrix
screen_mean = 0;
screen_std  = 1;
if (flag == 0)
    f           = screen_std*randn(n,n) + screen_mean;
    % enlarging the rough screen
    f           = repelem(f,n,n);
else
    if (flag == 1)
        f  = screen_std*rand(N,N);
    end
    
    if (flag ~= 1 && flag~= 0)
        error('Flag must be 0 or 1')
    end
end
% figure
% imagesc(f,'CDataMapping','scaled')
% title('2d roughness profile')

% expanding the roughness screen to match the gaussian by padding
% zeros if Nn < N
% cutting off the roughness screen if it exceeds the laser
if (flag == 0)
if(Nn ~= N)
    if (Nn < N)
        dN = N - Nn;
        f  = padarray(f,dN/2);
        f  = padarray(f',dN/2);
        disp('Screen smaller')
        
    else
        dN          = Nn - N;
        f(:,1:dN)   = [];
        f(:,end-dN+1:dN) = [];
        f(1:dN,:)   = [];
        f(end-dN+1:dN,:) = [];
        disp(2)
        
    end
end
end

% figure
% imagesc(f,'CDataMapping','scaled')
% title('2d roughness profile after moification')

% Multiplying beam with the rough screen
Psibeam = Psibeam.*exp(1i*kbeam.*f);


%% SLM generation with angular slits
% columns = length(Psibeam);
% rows    = length(Psibeam);
% alpha   = pi/6;% angular spacing between slits innermost
% phi     = pi/12; % width of each slit

mask    = make_angular_slit(Psibeam,d,alpha,phi);

kx  = kslm/10;
phi = atan(y./x);

T   = 0.5*(1-cos(kx*x - ll*phi));

% phase slm
T  = mod(ll*(atan2(y,x) + pi),2*pi);


T  = exp(1i*T);


%thresholding the slm
%T   = matrixthreshold(T);

% multiplying wiuth masl
T   = T.*mask;

% if (uu == 1)
%     figure
%     imagesc(T,'CDataMapping','scaled')
%     title('SLM')
%     %colormap(gray)
% end
%% Multiplying SLM with gaussian

G = T.*Psibeam;
% 
% figure
% imagesc(abs(G),'CDataMapping','scaled')
% title('Gaussian*SLM')

%fraunhaufer
[Uout,x2,y2] = fraunhofer_prop(G,lamb,d,Z);
MX = max(max((abs(Uout)).^2));
% figure
% imagesc(((abs(Uout)).^2)/MX,'CDataMapping','scaled')
% title(ll );



%% finding out the max rows and column value
[M,I] = max(abs(Uout(:)));
[I_row, I_col] = ind2sub(size(Uout),I);

% a 1d profile of the diffraction pattern
% figure
xx = x(I_row,:);
UU = abs(Uout(I_row,:)).^2;
% 
% plot(x(I_row,:),abs(Uout(I_row,:)));
% title(ll)

end



% % Angular Spectrum
% U  = Band_extended_angular_spect(T,x,y,d,lamb,Z);
% 
% figure
% imagesc(abs(U)/max(max(abs(Uout))),'CDataMapping','scaled')

























