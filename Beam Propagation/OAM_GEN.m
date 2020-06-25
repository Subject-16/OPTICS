clear all
close all

%code for generating OAM diffraction grating and OAM light
%35 mm plate and grating generated by interfering plane wave and beam with
%helical phase.All units in SI. Laser Beam not yet introduced

dscale = 1e-5;                %spacing in SI for generating x,y matrix 
scale  = -0.020:dscale:0.020;
I      = 0.66;%intensity of two waves kept same for simplicity
lamb   = 630*1e-9;% wavelength of 650nm
k      = 2*pi/lamb;
F      = 46.5; %Horizontal Adjustment factor
ll     = 1;%order of OAM
[x,y]  = meshgrid(scale);
Patt   = abs(I*exp(1i*k.*(x))+I*exp(1i*ll.*atan2(y,x))).^2;
Patt(3501,3501) = Patt(3501,3502);
figure;
image(Patt,'CDataMapping','scaled')
title('SLM grating with  L=1(to be refined further)')

M = Patt;
Mm = mean(Patt(:));
M(M>Mm)  = 1;
M(M~=1) = 0;
%image(M,'CDataMapping','scaled')
%generation of Diffraction pattern
Dz           = 20;%location of screen location from grating in m
[Uout,x2,y2] = fraunhofer_prop(M,lamb,dscale,Dz);
figure;

Inte = abs(Uout).^2;% intensity
image(log(Inte),'CDataMapping','scaled')% image plane
%plot(1:length(x2),Inte(y2==0))
title('Diffraction modes')
