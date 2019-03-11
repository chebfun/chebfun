function plot(f, varargin)
%PLOT Plot a BALLFUN on the ball
%
% PLOT(f, 'slices') plot a BALLFUN function and its slices on the planes 
% X-Y, Y-Z and X-Z
%
% See also QUIVER

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 1
    plot1(f)
elseif nargin == 2 && (strcmp(varargin{1},'slices') || strcmp(varargin{1},'slice'))
    plot2(f)
else
    error('CHEBFUN:BALLFUN:plot:input Invalid input arguments')
end
end

function plot1(f)
% PLOT Plot a BALLFUN function on the ballfun

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Define the size of F: 
[m,n,p] = size(f);

% m >= 25 and n, p >= 28
m = 25*(m < 25) + m*(m >= 25);
n = 28*(n < 28) + n*(n>=28);
p = 28*(p < 28) + p*(p>=28);

% Impose m = 1 [6] and n, p = 0 [4] to avoid errors in the plot
m = m + mod(1-mod(m,6),6);
n = n + mod(4-mod(n,4),4);
p = p + mod(4-mod(p,4),4);

% Get the coeffs of the ballfun function F(r,lam,th), lam is the
% azimuthal variable in [-pi,pi] and theta the polar variable in [0,pi]
F = coeffs3(f, m, n, p);

% Convert to values
ff = real(ballfun.coeffs2vals(F));

% Permute lambda and theta
ff = permute(ff,[1 3 2]);

% Evaluation points
r   = chebpts( m );
lam  = [pi*trigpts( n ); pi];
th = [pi*trigpts( p );pi]-pi/2;

% Remove doubled-up data
r = r(floor(m/2)+1:end);
th = th(floor(p/2)+1:end);

% Reverse theta : 1st element of the array is theta = pi (South Pole), last element is
% th = 0 (not included) (North Pole)
ff  = ff(floor(m/2)+1:end,[1 end:-1:floor(p/2)+1],:);
ff(:,:,end+1) = ff(:,:,1);

% Define the meshgrid
[tt, rr, ll] = meshgrid(th, r, lam);

% Slices in the cylinder to plot
% Find the indice of r = 0.5
[~,idr]=min(abs(r-0.5));
rslice = rr(idr,1,1);
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
hold off

axis([-1 1 -1 1 -1 1])
daspect([1 1 1])

camlight;
lighting phong;
material dull;

% Add label
xlabel('X')
ylabel('Y')
zlabel('Z')
end

function plot2(f)
% PLOT2 Plot a BALLFUN function on the ballfun and its slices

% Plot f on the plane X-Y
subplot(2,2,2);
plot(diskfun(f,'z'))
colorbar
xlabel('X')
ylabel('Y')

% Plot f on the plane X-Z
subplot(2,2,3);
plot(diskfun(f,'y'));
colorbar
xlabel('X')
ylabel('Z')

% Plot f on the plane Y-Z
subplot(2,2,4);
plot(diskfun(f,'x'));
colorbar
xlabel('Y')
ylabel('Z')

% Plot f
subplot(2,2,1);
plot1(f)
colorbar

set(gcf,'PaperPositionMode','auto','PaperPosition',[0 0 15 10])
end