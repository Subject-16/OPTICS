% replaces all Nan values with nearby values 
function [a] = replace_nan1(a)

nanIDX = find(isnan(a));
while(~isempty(nanIDX))
  a(nanIDX) = a(nanIDX+1);
  nanIDX      = find(isnan(a));
end

end
