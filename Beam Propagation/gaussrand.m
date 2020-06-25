%function for generating gaussian random variables as per requirement using
%the box mueller transform with zero mean nd unit variance
function y = gaussrand(n)%n is the number of gaussian rand nos required

w = 2;
y = zeros(1,n);
for ii = 1:floor(n/2)
  while w>=1
    x1 = 2*rand(1,1)-1;
    x2 = 2*rand(1,1)-1;
    w  = x1.^2+ x2.^2;
  end
  w  = sqrt((-2*log(w))/w);

  y1 = x1*w;
  y2 = x2*w;
  if(ii==1)
      y  = [y1,y2];
  else
      y  = [y,y1,y2];
      
  end
  
end


end
