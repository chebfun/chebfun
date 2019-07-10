function varargout = slice(f)
%SLICE   Plots slices (cross sections) of a BALLFUN.
%   SLICE(F) creates a slice plot of the BALLFUN at x = 0, y = 0 and z = 0.
%
% See also BALLFUN/PLOT.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check if the function is empty
if isempty(f)
    error('CHEBFUN:BALLFUN:slice:isempty','Function is empty.');
end

% Add a warning of the function is not real
if (f.isReal == 0)
    warning('CHEBFUN:BALLFUN:slice:isReal','Function is not real, plotting the real part.');
end

% Is the plot currently being held?
plotOnHold = ishold;
% Default plotting options
defaultOpts = {'facecolor', 'interp','edgecolor', 'none'};

% Define the size of f: 
[m,n,p] = size(f);

% m >= 25 and n, p >= 28
m = 25*(m < 25) + m*(m >= 25);
n = 28*(n < 28) + n*(n>=28);
p = 28*(p < 28) + p*(p>=28);

% Impose m = 1 [6] and n, p = 0 [4] to avoid errors in the plot
m = m + mod(1-mod(m,6),6);
n = n + mod(4-mod(n,4),4);
p = p + mod(4-mod(p,4),4);

% Discretization points
r = chebpts(m); r = r(floor(m/2)+1:end);
lam = linspace(-pi,pi,n);
th = linspace(0,pi,p);

% Evaluate the function at x = 0, y < 0
ff = permute(fevalm(f,r,-pi/2,th),[1 3 2]);
ff = real(ff);
% Plot the result
h = surf(zeros(length(r),length(th)),-r*sin(th),r*cos(th),ff,defaultOpts{:});

hold on

% Evaluate the function at x = 0, y > 0
ff = permute(fevalm(f,r,pi/2,th),[1 3 2]);
ff = real(ff);
% Plot the result
surf(zeros(length(r),length(th)),r*sin(th),r*cos(th),ff,defaultOpts{:})

% Evaluate the function at y = 0, x < 0
ff = permute(fevalm(f,r,pi,th),[1 3 2]);
ff = real(ff);
% Plot the result
surf(-r*sin(th),zeros(length(r),length(th)),r*cos(th),ff,defaultOpts{:})
% Evaluate the function at y = 0, x > 0
ff = permute(fevalm(f,r,0,th),[1 3 2]);
ff = real(ff);
% Plot the result
surf(r*sin(th),zeros(length(r),length(th)),r*cos(th),ff,defaultOpts{:})

% Evaluate the function at z = 0
ff = fevalm(f,r,lam,pi/2);
ff = real(ff);
% Plot the result
surf(r*cos(lam),r*sin(lam),zeros(length(r),length(lam)),ff,defaultOpts{:})

if ~plotOnHold
    hold off;
end

camlight('headlight');
lighting phong;
material dull;

axis([-1 1 -1 1 -1 1])
daspect([1 1 1])

if ( nargout > 0 )
    varargout = { h }; 
end
end