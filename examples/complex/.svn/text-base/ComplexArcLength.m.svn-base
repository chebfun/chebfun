%% Arc length in the complex plane
% Kuan Xu, 5th November 2012

%%
% (Chebfun example complex/ComplexArcLength.m)
% [Tags: #complex, #arclength, #roots]

%% 1. ARC LENGTH OF A CONTOUR
% Here is a contour discussed in one of the Example "A Keyhole
% Contour Integral" [1]: 

r = 0.2; R = 2; e = 0.1;
t = chebfun('t',[0 1]);                 % Parameter
c = [-R+e*1i -r+e*1i -r-e*1i -R-e*1i];
z = [ c(1) + t*(c(2)-c(1))              % Top of the keyhole
      c(2)*c(3).^t ./ c(2).^t           % Inner circle
      c(3) + t*(c(4)-c(3))              % Bottom of the keyhole
      c(4)*c(1).^t ./ c(4).^t];         % Outer circle
LW = 'LineWidth'; lw = 1.6;
plot(z,LW,lw), axis equal, xlim([-2.8 2.8])

%%
% The total length of the contour can be calculated with one command: 

L = arclength(z)

%%
% If the length of each patch is what you want to know, then you can make
% the input a quasimatrix whose columns are the different pieces.

z = [ c(1) + t*(c(2)-c(1))...           % Top of the keyhole
      c(2)*c(3).^t ./ c(2).^t...        % Inner circle
      c(3) + t*(c(4)-c(3))...           % Bottom of the keyhole
      c(4)*c(1).^t ./ c(4).^t];         % Outer circle

L = arclength(z)

%% 2. EQUIDISTRIBUTING POINTS ALONG A CONTOUR
% It is not uncommon that we need to discretize and sample over a 2D curve
% in the complex plane. For instance, we may want to uniformly distribute 
% points along the boundary of a domain when the boundary integral method is 
% used. That is, the points should be equally-spaced with respect
% to arc length. Consider the following star-shaped curve:

t = chebfun('t',[0 1]); 
s = exp(1i*2*pi*t).*(0.5*sin(8*pi*t).^2+0.5);
plot(s,LW,lw), axis equal, hold on

%%
% First, we calculate the total arc length of this closed curve.

L = arclength(s)

%%
% Suppose, for example, that we want to equidistribute 64 points.

N = 64;
h = L/N;
T = zeros(1,N);

%%
% We can locate the points by solving 63 rootfinding problems:
tic
len = cumsum(abs(diff(s)));
for k = 1:N-1
    T(k+1) = roots(len-k*h);    
end
toc

%%
% Now we have the coordinates of the points, and let's mark them on the
% curve.

P = s(T);
MS = 'MarkerSize'; ms = 20;
plot(P,'.r',MS,ms)

%%
% The rootfinding problems above, unfortunately, took a number of seconds,
% because the chebfun len is rather long:
length(len)

%%
% In principle, what we ought to do for this problem is first invert len 
% with the Chebfun commands INV or INV2, and then evaluate the result at 63
% equally spaced points.  However, INV and INV2 do not currently succeed at
% this, though INV2 can successfully invert on a subinterval, e.g. with the
% command g = inv2(len{0,0.1}). Further research is needed!

%%
% References:
%
% [1] Chebfun Example "A Keyhole Contour Integral" at
% http://www.chebfun.org/examples/complex/html/KeyholeContour.shtml