% code for making intensity vs OAM mode plot
% clc
% clear
close all
alph  = 0:0.2:pi/2;
del   = 0.25:0.25:4;
 



figure
% runinning the entire script for diff values of alpha
for al_ii = 1:length(del)

ll      = -10:1:10; % number of OAM modes needed
n       = 0;    % roughness index
thr     = 8.5e-2;
%thr     = ; %experimental threshold
alpha   = pi/12;% angular spacing between slits innermost
phi     = pi/12; % width of each slit

% 
% %plotting peaks
% for ii = 1:5:length(ll)
%     [x,Uout,xx,UU] = test1(n,ll(ii),ii,alpha,phi); % UU is abs(Uout).^2 and is 1d.
%     F2 = UU;
%     F1 = abs(Uout).^2;
% %     figure
% %     imagesc(F1,'CDataMapping','scaled')
% %     title('F1')
%     %close all % closing all plots
%     %% finding the peak values 
%     [Minim,MinIdx] = findpeaks(F2);
%     figure
%     plot(xx,F2);
%     title(ll(ii))
%     hold on
%     plot(xx(MinIdx),F2(MinIdx),'*') 
%     hold off
%     
% end


ll_value = zeros(1,length(ll));
for ii = 1:length(ll)
    [x,Uout,xx,UU] = test1(n,ll(ii),ii,alpha(1),phi,alph(4),del(al_ii));
    F2 = UU;
    %close all % closing all plots

    ll_value(ii) = F2(xx == 0); 
   
%     
end

% plotting the L dependence of OAM intensity
% figure
% plot(ll,ll_value,'-*')
% title('First second order Intensity')    


%% Obtaining the second order coherence function l and -l
LL      = ll;
thr2    = (thr^2)/2;% threshold 2nd order coherence peaks 
 

LL_value = zeros(1,length(LL));
for ii = 1:length(LL)
    
    [x1,U1,xx1,UU1] = test1(n,LL(ii),2,alpha,phi,alph(4),del(al_ii));
    [x2,U2,xx2,UU2] = test1(n,-LL(ii),2,alpha,phi,alph(4),del(al_ii));
    G2              = UU1.*UU2;
    
%     %%%%%%%%%%%%%%%%%%%% finding the peak values 
%     [Minim,MinIdx] = findpeaks(G2);
%     figure
%     plot(xx1,G2);
%     title(LL(ii))
%     hold on
%     plot(xx1(MinIdx),G2(MinIdx),'*')
%     hold off
%     %%%%%%%%%%%%%%%%%%%%%%%
    %close all % closing all plots

    LL_value(ii) = G2(xx1 == 0); 
    
%     
end



  

%% Section for making an analytical envelope 
%close all
clc
% First order correlation normalised



alpha1     = phi;
phi1       = alpha + phi;
Theo_corr1 = (alpha1/pi)^2*((sin(alpha1.*ll/2)./(alpha1.*ll/2)).^2)...
    .*(cos(phi1*ll/2)).^2;
facto      = max(ll_value)/max(Theo_corr1);
fshift     = max(ll_value) - max(Theo_corr1);

% figure
plot(ll,ll_value/max(ll_value),'-*')
title('First second order Intensity')    
hold on
plot(ll,Theo_corr1/max(Theo_corr1),'-')

% 
% % Second order correlation
% 
% clc
% 
% Theo_corr2 = (alpha1/pi)^4*((sin(alpha1.*ll/2)./(alpha1.*ll/2)).^4)...
%     .*(cos(phi1*ll*2))/(2*pi^6);
% 
% Theo_corr2 = Theo_corr1.^2;
% 
% figure
% plot(LL,LL_value/max(LL_value),'-*')
% title('Second second order Intensity')    
% hold on
% plot(LL,(Theo_corr2)/max(Theo_corr2),'-')
% hold off

end

%% Theoretical for OV beam + guassian SIM
alpha   = pi/12;% angular spacing between slits innermost
phi     = pi/12; % width of each slit

ll    = -10:0.5:10; 
alph  = 0:0.2:1.6;
del   = 0.25:0.25:4;
w     = .35e-3;
w1    = sqrt(2)*(w);
chg   = 1;
a     = sin(alph(4));
b     = cos(alph(4));
R1    = 1/(2*pi)*(a^2*w^2)/4;
R2    = 1/(2*pi)*(b^2*w1^2)/4;
R3    = 1/(2*pi)*(a*b*w^2*w1^2)/(2*(w1^2 + w^2));
al_ii = 4;
a_alph=  4;

figure
for al_ii = 1:length(del)
    
alpha1  = phi;
phi1    = alpha + phi;
% 
% a     = sin(alph(a_alph));
% b     = cos(alph(a_alph));
% R1    = 1/(2*pi)*(a^4*w^4)/4;
% R2    = 1/(2*pi)*(b^2*w1^2)/4;
% R3    = 1/(2*pi)*(a*b*w^2*w1^2)/(2*(w1^2 + w^2));

Theo       = R1*2./(ll.^2).*((sin(ll*alpha/2)).^2).*(1 + cos(ll*phi1)) + R2*2./((ll + chg).^2).*((sin((ll + chg)*alpha/2)).^2).*(1 + cos((ll + chg)*phi1))...
    + R3*2./((ll + chg).*ll).*sin(ll*alpha/2).*sin((ll + chg)*alpha/2).*(cos(del(al_ii)) + cos(ll*phi1 - del(al_ii)) + cos(chg*phi1 + del(al_ii))...
    + cos((ll + chg)*phi1 + del(al_ii)) );

Theo1      = R1*2./(-ll.^2).*((sin(-ll*alpha/2)).^2).*(1 + cos(-ll*phi1)) + R2*2./((-ll + chg).^2).*((sin((-ll + chg)*alpha/2)).^2).*(1 + cos((-ll + chg)*phi1))...
    + R3*2./((-ll + chg).*-ll).*sin(-ll*alpha/2).*sin((-ll + chg)*alpha/2).*(cos(del(al_ii)) + cos(-ll*phi1 - del(al_ii)) + cos(chg*phi1 + del(al_ii))...
    + cos((-ll + chg)*phi1 + del(al_ii)) );

Theo_2nd_order = Theo1.*Theo;

hold on
plot(ll,Theo_2nd_order/max(abs(Theo_2nd_order)),'*-')
%plot(ll,Theo/max(abs(Theo)),'*-')


end















% T  = mod(ll*(atan2(y,x) + pi),2*pi);
% T  = exp(1i*T);








    
    