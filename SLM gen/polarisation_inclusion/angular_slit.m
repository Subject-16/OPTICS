
%% Make pie sector.
xCenter = columns/2;
yCenter = rows/2;
theta1 = (pi/4-pi/24) : 0.001 : (pi/4+pi/24);
theta2 = theta1 + pi;
radius = min([rows/2, columns/2]);

x1 = radius * cos(theta1) + xCenter;
y1 = radius * sin(theta1) + yCenter;
x2 = radius * cos(theta2) + xCenter;
y2 = radius * sin(theta2) + yCenter;

x  = [x2]; y = [y2];

% Add point at center to beginning and end of array and shift to center.
x = [xCenter, x, xCenter];
y = [yCenter, y, yCenter];
subplot(2, 2, 2);
plot(x, y);
title('Sector Plot', 'FontSize', fontSize);
axis square;
grid on;

mask = poly2mask(x,y, size(rgbImage, 1), size(rgbImage, 2));
% Display the mask image.
subplot(2, 2, 3);
imshow(mask);
title('Mask Image', 'FontSize', fontSize);

%% Mask the image using bsxfun() function
maskedRgbImage = bsxfun(@times, rgbImage, cast(mask, class(rgbImage)));
% Display the mask image.
subplot(2, 2, 4);
imshow(maskedRgbImage);
title('Masked Color Image', 'FontSize', fontSize);