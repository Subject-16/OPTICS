%% SCript for OAM mode analysis

function [MM,l_weight] = OAM_decomp_LG_basis(alph1,Vin,w1,x,y,d,Z,lamb,g)

% Need to create a basis of l states from -1 to +1 
MM       = -10:1:10;
l_weight = zeros(size(MM));
r        = sqrt(x.^2 + y.^2);

%V_areastore = zeros(size(MM));

for ii = 1:length(MM)
% 
% TT    = mod(MM(ii)*(atan2(y,x) + pi),2*pi);
% T     = exp(1i*TT);
% Vbeam = exp(-(x.^2 + y.^2)/w1^2).*T;
% 
% [V_eigen,~,~] = fraunhofer_prop(Vbeam,lamb,d,Z); % creating OV beam with desired l value
% 
% MX = max(max((abs(V_eigen)).^2));

%% Forming the LG basis(p = 0)

TT    = mod(MM(ii)*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);
p     = 0;

Vbeam = sqrt(2*factorial(p)./(pi*w1*factorial(abs(MM(ii) + p)))).*(r.*sqrt(2)./w1).^(abs(MM(ii)))...
    .*exp(-((r.^2)./w1^2)).*T;
Vbeam = exp(-((r.^2)./w1^2)).*T;

%% Normalizing LG beam


%V_areastore(ii) = V_area;

V_eigen =  (Vbeam)/max(abs((Vbeam(:))));
V_area  = trapz(g,trapz(g,(abs(V_eigen)).^2,2));
%V_eigen = V_eigen;

% figure
% imagesc(((abs(V_eigen)).^2),'CDataMapping','scaled')
% title(MM(ii))
% ylabel(V_area)

%% Propagating the beam and finding the central field value

U           = conj(V_eigen).*Vin; % performing the inner product
[U_out,~,~] = fraunhofer_prop(U,lamb,d,Z); % propagating the resulting beam
% 
% figure
% imagesc(((abs(U_out)).^2),'CDataMapping','scaled')
%  

m = (x == 0).*(y == 0);
m = logical(m); % finding logically central index

l_weight(ii) = abs(U_out(m)).^2; % output is the central value of far field

end
l_weight = l_weight/(sum(l_weight));
%
figure
bar(MM,l_weight);
title(alph1)

% need to plot l value corresponding to peak vs alpha/del
[Maxima_val,MaximaIdx] = findpeaks(l_weight); 

g1 = 2;
[a,b] = max(Maxima_val);

b = MaximaIdx(b);
a = MM(b);
end

