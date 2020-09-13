% code for making intensity vs OAM mode plot
clc
clear
%close all

ll      = -15:0.2:15; % number of OAM modes needed
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
    [x,Uout,xx,UU] = test1(n,ll(ii),ii,alpha,phi);
    F2 = UU;
    %close all % closing all plots
    %% finding the peak values and averaging
     [Minim,MinIdx] = findpeaks(F2);

    a   = xx(MinIdx);
    den = 0;
    mean_OAM = 0; % will contain the mean of the peaks
    %% finding the peaks corresponding to the OAM mode
    for jj = 1:length(a)
        if (abs(a(jj))== 0) %&& abs(a(jj))>0.5e-3)
            %if( Minim(jj) > thr) % setting cutoff for peak selection
                mean_OAM = Minim(jj) + mean_OAM;
                den      = den + 1;
                
            %end
        end
        
    end
    mean_OAM     = mean_OAM/den;
    ll_value(ii) = mean_OAM;
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
    
    [x1,U1,xx1,UU1] = test1(n,LL(ii),2,alpha,phi);
    [x2,U2,xx2,UU2] = test1(n,-LL(ii),2,alpha,phi);
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
    % finding the peak values and averaging
     [Minim,MinIdx] = findpeaks(G2);

    a   = xx1(MinIdx);
    den = 0;
    mean_OAM = 0; % will contain the mean of the peaks
    % finding the peaks corresponding to the OAM correlation mode
    for jj = 1:length(a)
        if (abs(a(jj)) == 0) % && abs(a(jj))>0.49e-3)
            %if( Minim(jj) > thr2) % setting cutoff for peak selection
                mean_OAM = Minim(jj) + mean_OAM;
                den      = den + 1;
                
           % end
        end
        
    end
    mean_OAM     = mean_OAM/den;
    LL_value(ii) = mean_OAM; 
    
%     
end

% % plotting the L dependence of OAM intensity
% figure
% plot(LL,LL_value,'-*')
% title('Second order coherence function')   



% %% Binning values together to only leave integer values
% 
% 
% % multiplying the intensities(E^2) (ensemble average not being done as ensembles
% % are all coherent
% 
% bin       = ll(1):1:ll(end);
% o         = discretize(ll,bin);
% ll_binned = zeros(1,length(bin));% new vector for summed values
% LL_binned = zeros(1,length(bin));
% 
% for ii = 1:length(o)
%     ll_binned(o(ii)) = ll_value(ii) + ll_binned(o(ii));
%     LL_binned(o(ii)) = LL_value(ii) + LL_binned(o(ii));
% end
% 
% 
% figure
% plot(bin,ll_binned,'-*')
% title('First order coherence function Binned')   
% 
% figure
% plot(bin,LL_binned,'-*')
% title('Second order coherence function Binned')   

%% Section for making an analytical envelope 
close all
clc
% First order correlation normalised



alpha1     = phi;
phi1       = alpha + phi;
Theo_corr1 = (alpha1/pi)^2*((sin(alpha1.*ll/2)./(alpha1.*ll/2)).^2)...
    .*(cos(phi1*ll/2)).^2;
facto      = max(ll_value)/max(Theo_corr1);
fshift     = max(ll_value) - max(Theo_corr1);

figure
plot(ll,ll_value/max(ll_value),'-*')
title('First second order Intensity')    
hold
plot(ll,Theo_corr1/max(Theo_corr1),'-')


% Second order correlation

clc

Theo_corr2 = (alpha1/pi)^4*((sin(alpha1.*ll/2)./(alpha1.*ll/2)).^4)...
    .*(cos(phi1*ll*2))/(2*pi^6);

Theo_corr2 = Theo_corr1.^2;

figure
plot(LL,LL_value/max(LL_value),'-*')
title('Second second order Intensity')    
hold
plot(LL,(Theo_corr2)/max(Theo_corr2),'-')



















% T  = mod(ll*(atan2(y,x) + pi),2*pi);
% T  = exp(1i*T);








    
    