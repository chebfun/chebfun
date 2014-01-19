%% Implicit functions and fun with contours
% Stefan Güttel, 30th July 2012

%%
% (Chebfun example fun/ContourFun.m)
% [Tags: #implicitfunction, #fun, #contour]

close all
clear all
LW = 'linewidth'; lw = 2.5;
MS = 'markersize'; ms = 10;

%%
% Matlab provides a function called |contours|, which returns a vector of
% points on level lines formed by the entries of a matrix. This function 
% can easily be used to define chebfuns satisfying an implicit equation, at 
% least approximately. To demonstrate this, suppose we want to create a 
% chebfun that corresponds to a circle of radius 0.8 centered at the 
% origin. The points on this circle satisfy z(x,y) = 0 with 

z = @(x,y) x.^2 + y.^2 - 0.8.^2;

%%
% We now sample z at 20 x 20 grid points on the square [-1,1] x [1,1], and 
% use Matlab's |contours| to detect the level lines of the circle:

x = linspace(-1,1,20);
[X,Y] = meshgrid(x,x);
Z = z(X,Y);
C = contours(X,Y,Z,[0,0]);

%%
% The following plot shows that the algorithm does quite a good job by
% actually interpolating between the grid points:

plot(.8*exp(2i*pi*linspace(0,1,200)),'k-',LW,lw)
hold on, axis square
plot(X(:),Y(:),'b+',MS,ms)
plot(C(1,2:end),C(2,2:end),'ro',MS,ms)


%%
% The following lines will convert the coordinates of the contours to a 
% quasimatrix f which each column corresponding to a contour (in this case,
% there is only one contour):

j = 1; f = chebfun();
while j < length(C),
    k = j + C(2,j); D = C(:,j+1:k); j = k + 1;
    f = [ f , chebfun(D(1,:)+1i*D(2,:)) ]; 
end
figure
plot(f,'r',LW,lw)
axis square, grid on

%%
% Of course, we cannot expect this to be a perfect circle as we have only
% sampled at 20 grid points in each coordinate direction, but for some
% practical purposes this accuracy may be sufficient:

norm(z(real(f),imag(f)),inf)

%%
% One example where high precision is certainly not necessary is when a
% function has do be defined from an image. To demonstrate this, let's read 
% an image showing a whiteboard in our visitors office. I'm grateful to 
% Yuji Nakatsukasa for taking this snapshot. 

A = imread('stefun.jpg');
imshow(A)
axis on, hold on

%%
% In the following we will actually operate on a grayscale version of this
% image, which is computed by the following two lines:

A = 0.299*double(A(:,:,1)) + 0.587*double(A(:,:,2)) + 0.114*double(A(:,:,3));
A = A/max(max(A));

%%
% We crop the image to the region of the graph on the upper right, and
% compute contour lines using the |contours| command:

x = 185:463;
y = 60:210;
B = A(y,x);
[X,Y] = meshgrid(x,y);
C = contours(X,Y,B,[0.8,0.8]);

%%
% As above, the following lines will convert the coordinates of the 
% contours to a quasimatrix f:

j = 1; f = chebfun();
while j < length(C),
    k = j + C(2,j); D = C(:,j+1:k); j = k + 1;
    f = [ f , chebfun(D(1,:)+1i*D(2,:)) ]; 
end

%%
% It turns out that f has four columns, and there are pairs of chebfuns 
% with function values very close to each other. These pairs correspond to 
% the upper and lower edge of a single line on the whiteboard. 
% By taking the arithmetic mean of two neighbouring chebfuns we obtain a 
% pretty good agreement with the curves on the board:

xaxis = (f(:,1) + f(:,2))/2;
curve = (f(:,3) + f(:,4))/2;
plot([xaxis,curve],LW,lw)

%%
% Since these are chebfuns, we can do all kinds of things with them. For
% example, find the area enclosed by the x-axis and the curve. Note that 
% the real part  of each chebfun will correspond to the x-coordinate, and 
% the imaginary part is the y-coordinate measured in pixels. The unit of 
% area therefore is square-pixels.

meanaxis = 1i*mean(imag(xaxis));
curve2 = curve - meanaxis;
Area = abs(sum(imag(curve2).*real(curve2)/2));
sprintf('The enclosed area is %0.1f square-pixels.',Area)

%%
% Or, we could compute a degree 3 best polynomial approximation to the 
% curve, overlay it with the plot, and show the intersections of both
% functions: 

p = 1i*remez(imag(curve),3) + real(curve);
plot(p,'r--',LW,lw)
rts = roots(curve - p);
plot(p(rts),'ro',MS,ms)

%%
% One may wonder, how close is the lower curve on the whiteboard to a 
% circle? In order to answer this question we again crop to the region of
% interest, and compute contours, which in this case gives a quasimatrix 
% with two columns corresponding to the edges of the painted curve. Let's 
% take the first column only as a prediction to the curve:

x = 310:470;
y = 240:410;
B = A(y,x);
[X,Y] = meshgrid(x,y);
C = contours(X,Y,B,[0.8,0.8]);
j = 1; f = chebfun();
while j < length(C),
    k = j + C(2,j); D = C(:,j+1:k); j = k + 1;
    f = [ f , chebfun(D(1,:)+1i*D(2,:)) ]; 
end
curve = f(:,1);
plot(curve,'m',LW,lw)

%%
% We can compute the area and centroid of this curve (see also [1]), and 
% draw a perfect circle with the same area for comparison:

Area = sum(real(curve).*diff(imag(curve)));
centroid = sum(diff(curve).*curve.*conj(curve))/(2i*Area);
radius = sqrt(Area/pi);
plot(centroid,'m+',MS,ms)
circle = chebfun('exp(1i*pi*x)');
circle = radius*circle + centroid;
plot(circle,'k--',LW,lw)

%%
% Last and least, we go crazy and approximate ourselfs via a collection of
% quasimatrices plotted with three shades of gray:

x = 10:170;
y = 90:310;
B = A(y,x);
[X,Y] = meshgrid(x,y);
levels = [.2 .4 .8];
C = contours(X,Y,B,levels);
j = 1; f = chebfun(); level = [];
while j < length(C),
    k = j + C(2,j);  D = C(:,j+1:k); level = [ level , C(1,j) ];
    j = k + 1; f = [ f , chebfun(D(1,:)+1i*D(2,:)) ];
end
figure
for j = 1:length(levels),
    ind = find(level == levels(j));
    plot(f(:,ind),'Color',[levels(j) levels(j) levels(j)],LW,lw)
    hold on
end
axis ij, axis image

%%
% I refrain from doing any calculations on this!

%%
% References:
%
% [1] http://www2.maths.ox.ac.uk/chebfun/examples/geom/html/Area.shtml
%


