% Fucntion for making a sector mask for angular slits
function [mask] = make_angular_slit(Psibeam,d,alpha,phi) % beam input for length and d for spacing in theta
%% Make pie sector.
columns = length(Psibeam);
rows    = length(Psibeam);
% alpha   = pi/6;% angular spacing between slits innermost
% phi     = pi/12; % width of each slit

xCenter  = columns/2;
yCenter  = rows/2;
starang1 = pi/2-(alpha/2 + phi); 
theta    = starang1: d: starang1 + phi ;
theta1   = theta + pi;
theta2   = theta + (alpha + phi);
theta2   = theta2 + pi;
radius   =  min([rows/2, columns/2]);

x1 = radius * cos(theta1) + xCenter;
y1 = radius * sin(theta1) + yCenter;
x2 = radius * cos(theta2) + xCenter;
y2 = radius * sin(theta2) + yCenter;

%% Making first slit
x  = [x1]; y = [y1];

% Add point at center to beginning and end of array and shift to center.
x = [xCenter, x, xCenter];
y = [yCenter, y, yCenter];

mask = poly2mask(x,y, rows, columns);
% 
%% Making second slit
x  = [x2]; y = [y2];

% Add point at center to beginning and end of array and shift to center.
x = [xCenter, x, xCenter];
y = [yCenter, y, yCenter];

mask = mask + poly2mask(x,y, rows, columns);



% %% Making a single slit comment this section 
% 
% starang1 = pi/2 - (phi/2);
% starang2 = pi/2 + (phi/2);
% 
% theta    = starang1: d: starang2 ;
% theta    = theta + pi;
% x1       = radius * cos(theta) + xCenter;
% y1       = radius * sin(theta) + yCenter;
% 
% x  = [x1]; y = [y1];
% 
% % Add point at center to beginning and end of array and shift to center.
% x = [xCenter, x, xCenter];
% y = [yCenter, y, yCenter];
% 
% mask = poly2mask(x,y, rows, columns);
%Display the mask image.
% figure
% imagesc(mask,'CDataMapping','scaled');
% title('Mask Image');
