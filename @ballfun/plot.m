function plot(f,varargin)
% PLOT Plot a BALLFUN function on the ballfun

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the coeffs of the ballfun function F(r,lam,th), lam is the
% azimuthal variable in [-pi,pi] and theta the polar variable in [0,pi]
F = f.coeffs;

% Expand the tensor if the number of coeffs in r, lam or th isn't equal to
% 0 modulo 4 to avoid errors in the plot
% or if the number of coefficients is less than 20
G = expand(F);
[m,n,p] = size(G);

% Permute lambda and theta
G = permute(G,[1 3 2]);

ff = real(ballfun.coeffs2vals(G));

r   = chebpts( m );
lam  = [pi*trigpts( n ); pi];
th = pi*trigpts( p )-pi/2;

% Remove doubled-up data
r   = r(floor(m/2)+1:end);
th = th(floor(p/2)+1:end);
% Reverse theta : 1st element of the array is theta = pi (South Pole), last element is
% th = 0 (not included) (North Pole)
ff  = ff(floor(m/2)+1:end,[1 end:-1:floor(p/2)+2],:);
ff(:,:,end+1) = ff(:,:,1);

[tt, rr, ll] = meshgrid(th, r, lam);

% Slices in the cylinder to plot
rslice = rr(floor(m/4),1,1);
tslice = tt(1,[1,floor(p/4)+1],1);
lslice = ll(1,1,[1,floor(n/4)+1]);

hslicer = slice(tt,rr,ll,ff,tslice,rslice,lslice);
hold on
for j = 1:numel(hslicer)
    h = hslicer(j);
    [xs,ys,zs] = sph2cart(h.ZData,h.XData,h.YData);
    surf(xs,ys,zs,h.CData,'EdgeColor','none','FaceColor','Interp');
end
delete(hslicer);
axis([-1 1 -1 1 -1 1])
daspect([1 1 1])

% Sphere grid
if (nargin==1 || varargin{1}~="Grid")
    [x, y, z] = sphere(16);
    surf(x, y, z, 'linestyle',':','FaceAlpha', 0);
end
hold off

camlight;
lighting phong;
colormap(jet);
colorbar;

% Add label
xlabel('X')
ylabel('Y')
zlabel('Z')
end

function G = expand(F)

S = size(F);
% Avoid an issue if p = 1
if length(S) == 2
    S(3) = 1;
end
m = S(1); n = S(2); p = S(3);

% Compute the new size of the tensor
mExpand = max(m + mod(4-mod(m,4),4), 20);
nExpand = max(n + mod(4-mod(n,4),4), 20);
pExpand = max(p + mod(4-mod(p,4),4), 20);

% Find the list of coefficients which corresponds to the coefficients of F
if mod(n,2) == 0
    ListLambda = 1+floor((nExpand-n)/2):n+floor((nExpand-n)/2);    
else
    ListLambda = 1+ceil((nExpand-n)/2):n+ceil((nExpand-n)/2); 
end

if mod(p,2) == 0
    ListTheta = 1+floor((pExpand-p)/2):p+floor((pExpand-p)/2);    
else
    ListTheta = 1+ceil((pExpand-p)/2):p+ceil((pExpand-p)/2); 
end

% Return the expansion G of F
G = zeros(mExpand, nExpand, pExpand);
G(1:m, ListLambda, ListTheta) = F;
end
