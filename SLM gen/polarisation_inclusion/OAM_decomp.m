%% SCript for OAM mode analysis

function [MM,l_weight] = OAM_decomp(Vin,w1,x,y,d,Z,lamb,g)

% Need to create a basis of l states from -1 to +1 
MM       = -5:0.1:5;
l_weight = zeros(size(MM));

for ii = 1:length(MM)

TT    = mod(MM(ii)*(atan2(y,x) + pi),2*pi);
T     = exp(1i*TT);
Vbeam = exp(-(x.^2 + y.^2)/w1^2).*T;

[V_eigen,~,~] = fraunhofer_prop(Vbeam,lamb,d,Z);

MX = max(max((abs(V_eigen)).^2));

% figure
% imagesc(((abs(V_eigen)).^2)/MX,'CDataMapping','scaled')

% V_area  = trapz(g,trapz(g,abs(V_eigen),2));
% 
% V_eigen =  (V_eigen)/V_area;

U = conj(V_eigen).*Vin;
[U_out,~,~] = fraunhofer_prop(U,lamb,d,Z);

% figure
% imagesc(((abs(U_out)).^2)/MX,'CDataMapping','scaled')


m = (x == 0).*(y == 0);
m = logical(m); % finding logically central index

l_weight(ii) = abs(U_out(m));

end

figure
bar(MM,l_weight);

end

