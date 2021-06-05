% code for making intensity vs OAM mode plot
% clc
% clear
clear
clc
close all
alph  = 0:0.2:pi/2;
del   = 0:1:4;
 
LLK = 1;% del or alpha index when other is varied

ll      = -20:0.1:20; % number of OAM modes needed
ll(floor(length(ll)/2)+1) = 0.002;

M_val        = -10:1:10;
ll_value_mat = zeros(length(alph),length(ll));
OAM_decomp   = zeros(length(alph),length(M_val));

figure
% running the entire script for diff values of alpha
for al_ii = 1:length(alph)

    
    
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
    if(ii == 1)
    [M_val,OAM_decomp(al_ii,:),x,Uout,xx,UU] = test1(n,ll(ii),ii,alpha,phi,alph(al_ii),del(LLK));
    
    else
    [~,~,x,Uout,xx,UU] = test1(n,ll(ii),ii,alpha,phi,alph(al_ii),del(LLK));
    end
    F2 = UU;
    %close all % closing all plots

    ll_value(ii) = F2(xx == 0); 
   
%     
end

ll_value_mat(al_ii,:) = ll_value;


% plotting the L dependence of OAM intensity
% figure
% plot(ll,ll_value,'-*')
% title('First second order Intensity')    


%% Obtaining the second order coherence function l and -l
LL      = ll;
thr2    = (thr^2)/2;% threshold 2nd order coherence peaks 
 

LL_value = zeros(1,length(LL));
for ii = 1:length(LL)
    
    [~,~,x1,U1,xx1,UU1] = test1(n,LL(ii),2,alpha,phi,alph(al_ii),del(LLK));
    [~,~,x2,U2,xx2,UU2] = test1(n,-LL(ii),2,alpha,phi,alph(al_ii),del(LLK));
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

figure
plot(ll,LL_value)
  


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


alph  = 0:0.2:pi/2;
del   = 0:1:4;
LLK = 1;% del or alpha index when other is varied
ll      = -20:0.1:20; % number of OAM modes needed
ll(floor(length(ll)/2)+1) = 0.002;



alpha   = pi/12;% angular spacing between slits innermost
phi     = pi/12; % width of each slit
ll_Theo_value_mat = zeros(length(alph),length(ll));

%% Section for making an analytical envelope 
% First order correlation normalised




alpha1     = phi;
phi1       = alpha + phi;
Theo_corr1 = (alpha1/pi)^2*((sin(alpha1.*ll/2)./(alpha1.*ll/2)).^2)...
    .*(cos(phi1*ll/2)).^2;

% % figure
% plot(ll,ll_value/max(ll_value),'-*')
% title('First second order Intensity')    
% hold on
% plot(ll,Theo_corr1/max(Theo_corr1),'-')

% 

%alph  = 0:0.2:1.6; defined at the beggining
%del   = 0.25:0.25:4;


w     = 1; % no need for exact values for theoretical plot
w1    = (w);
chg   = -1.0001; % charge of OV beam
% a     = sin(alph(4));
% b     = cos(alph(4));
% R1    = 1/(2*pi)*(a^2*w^2)/4;
% R2    = 1/(2*pi)*(b^2*w1^2)/4;
% R3    = 1/(2*pi)*(a*b*w^2*w1^2)/(2*(w1^2 + w^2));
% al_ii = 4;
% a_alph=  4;

%figure
for al_ii = 1:length(alph)
    
 
alpha1  = phi;
phi1    = alpha + phi;
% 

a_alph= al_ii;
a     = sin(alph(a_alph));
b     = cos(alph(a_alph));
R1    = 1/(2*pi)*(a^4*w^4)/4;
R2    = 1/(2*pi)*(b^2*w1^2)/4;
R3    = 1/(2*pi)*(a*b*w^2*w1^2)/(2*(w1^2 + w^2));

Theo       = R1*2./(ll.^2).*((sin(ll*alpha/2)).^2).*(1 + cos(ll*phi1)) + R2*2./((ll + chg).^2).*((sin((ll + chg)*alpha/2)).^2).*(1 + cos((ll + chg)*phi1))...
    + R3*2./((ll + chg).*ll).*sin(ll*alpha/2).*sin((ll + chg)*alpha/2).*(cos(del(LLK)) + cos(ll*phi1 - del(LLK)) + cos(chg*phi1 + del(LLK))...
    + cos((ll + chg)*phi1 + del(LLK)) );

Theo1      = R1*2./(-ll.^2).*((sin(-ll*alpha/2)).^2).*(1 + cos(-ll*phi1)) + R2*2./((-ll + chg).^2).*((sin((-ll + chg)*alpha/2)).^2).*(1 + cos((-ll + chg)*phi1))...
    + R3*2./((-ll + chg).*-ll).*sin(-ll*alpha/2).*sin((-ll + chg)*alpha/2).*(cos(del(LLK)) + cos(-ll*phi1 - del(LLK)) + cos(chg*phi1 + del(LLK))...
    + cos((-ll + chg)*phi1 + del(LLK)) );

Theo_2nd_order = Theo1.*Theo;

% hold on
% %plot(ll,Theo_2nd_order/max(abs(Theo_2nd_order)),'*-')
% plot(ll,abs(Theo)/max(abs(Theo)),'*-')
% plot(ll,ll_value/max(ll_value),'-*')

ll_Theo_value_mat(al_ii,:) = Theo;

end


% 
% 
% % 
% figure for showing shifts in OAM plots with varying alpha
figure
for ii = 1:length(alph)
    subplot(2,4,ii)
    plot(ll,ll_Theo_value_mat(ii,:)/max(ll_Theo_value_mat(ii,:)))
    hold on
    plot(ll,Theo_corr1/max(Theo_corr1))
    xlabel('l value')
    ylabel('OAM intensity')
    xlim([-20,20])
    ylim([0,1])
    x = alph(ii)*180/pi;
    title(sprintf('%0.1f', x))
    grid on
    hold off
end
suptitle('Shift in the OAM intensity(theoretical) with change in \alpha')

% 
% 
% 
% %% To deduce the lateral shift from the gaussian coherence
figure
plot(ll,Theo_corr1/max(Theo_corr1))
hold on
%plot(ll,(Theo)/max(abs(Theo)))
plot(ll,ll_value_mat(2,:)/max(ll_value_mat(2,:)))
plot(ll,ll_Theo_value_mat(2,:)/max(ll_Theo_value_mat(2,:)))
xlim ([-20,20])
ylim([0,1])
hold off


l_shift_theo = zeros(1,length(alph));
l_shift_exp  = zeros(1,length(alph));
for ii = 1:length(alph)
    
[a,b] = max(ll_Theo_value_mat(ii,:));% a is maximum value. b is the index
[A,B] = max(ll_value_mat(ii,:));% A is maximum value. B is the index
l_shift_theo(ii) = ll(b);
l_shift_exp(ii)  = ll(B);

end

decomp_ratio = zeros(size(alph));

% Modal decomp vs alpha
figure
for ii = 1:length(alph)
    subplot(2,4,ii)
    bar(M_val,OAM_decomp(ii,:))
    x = alph(ii)*180/pi;
    title(sprintf('%0.1f', x))
    xlabel('OAM value')
    xlim([-5,5])
    ylabel('Weights')
    
    % Find the ratio of the 0 and 1 components
    OA = OAM_decomp(ii,:);
    decomp_ratio(ii) = OA(M_val == -1);
end
suptitle('Modal Decomposition vs Input Polarisation')


% exper_theo shifts vs alpha
figure
plot(alph*180/pi,l_shift_theo,'-*')
hold on
plot(alph*180/pi,l_shift_exp,'-o')
hold on
plot(alph*180/pi,decomp_ratio,'-x')
xlabel('\alpha of the input polarisation')
ylabel('Shifts of the Theoretical and Simulation OAM Intensities ')



% Second order plots for OV beam for alpha = 80.2 degrees(last)
figure
plot(ll,LL_value)
xlabel('l value')
ylabel('Correlation Value(I_{l - l})')
grid on



% T  = mod(ll*(atan2(y,x) + pi),2*pi);
% T  = exp(1i*T);








    
    