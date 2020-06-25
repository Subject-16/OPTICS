function z = rectslit(x,y,d)% d is the width of slit
   z   = zeros(size(x));
   mid = floor(length(z));% position of mid point
   
   z(abs(x)<d/2) = 1;
   z(x==mid+d/2)=0.5;
   z(x==mid-d/2)=0.5;
end

   
