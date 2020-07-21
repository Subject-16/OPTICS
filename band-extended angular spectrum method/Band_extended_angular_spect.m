% Extended Band Angular Spectrum Method Padded input required
% for the 2d case first 1d ffft of the columns are performed. Then the rows
% then the shifting is done.

function [t_3_2] =    Band_extended_angular_spect(t,x,y,Dx,lamb,Z)

n       = length(x);
fc      = 1/(lamb)*1/sqrt(1 + (Z/(n*Dx))^2); % theoretical fourier limit on local fourier frequency
k       = (-n/2:n/2-1)/(n*Dx);
[kx,ky] = meshgrid(k,k);


iflag = -1;                 %for NUFFT
eps   = 10^(-12);           % accuracy of NUFFT
K     = n/2/max(max(kx));
fcn   = 1/2*sqrt(n/lamb/Z); % f_extend
ss    = fcn/max(max(abs(kx)));
zc    = n*Dx^2/(lamb)*sqrt(1-(lamb/(2*Dx))^2)
if Z < zc
    kxn = kx;
    kyn = ky;
else
    kxn = kx*(ss-0.0);
    kyn = ky*(ss-0.0);
end

p  = sqrt(1/lamb^2- kxn.^2 - kyn.^2);
G  = exp(1i*2*pi*Z.*p);


disp('band-extended ASM:')
tic
%% doing nufft step by step
%preallocations
t_asmNUFT1 = zeros(size(t));
t_asmNUFT2 = zeros(size(t));

t_3_1      = zeros(size(t));
t_3_2      = zeros(size(t));

max_x     = max(max(abs(x)));

%% NON-UNIFORM FOURIER TRANSFORM
%columns
for ii = 1:length(x)
   t_asmNUFT1(:,ii) = nufft1d3(n,y(:,ii)/max_x*pi,t(:,ii),iflag,eps,n,kyn(:,ii)*K);
end
%rows
for ii = 1:length(x)
   t_asmNUFT2(ii,:) = nufft1d3(n,x(ii,:)/max_x*pi,t_asmNUFT1(ii,:),iflag,eps,n,kxn(ii,:)*K);
end
%% multiplying phase factor
G  = G.*t_asmNUFT2;



%% NON-UNIFORM INVERSE FOURIER TRANSFORM
%columns
for ii = 1:length(x)
    t_3_1(:,ii) = nufft1d3(n,y(:,ii)/max_x*pi,G(:,ii),-iflag,eps,n,kyn(:,ii)*K);
end
%rows
for ii = 1:length(x)
    t_3_2(ii,:) = nufft1d3(n,x(ii,:)/max_x*pi,t_3_1(ii,:),-iflag,eps,n,kxn(ii,:)*K);
end

toc



%%
t_3_2 = t_3_2/max(max((abs(t_3_2))));
% t_3 = t_3(n/2-n/4+1:n/2+n/4,1);
% phase_asm_ex = (angle(t_3));
% amplitude_asm_ex = abs(t_3);
% 


%figure,plot(x(n/2-n/4+1:n/2+n/4),(amplitude_asm_ex));title('band-extended ASM amplitude ')
% figure,plot(x(n/2-n/4+1:n/2+n/4),(phase_asm_ex));title('band-extended ASM phase ')


end