function [pol] = laguerre_LG(l,p,r)

  pol = 0; 
  for m = 0:p
      pol = pol + ((-1)^m)*nchoosek(p+l,p-m)*(r.^m)./factorial(m);
    
  end
  
  
end