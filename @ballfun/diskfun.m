function g = diskfun(f, varargin)
% DISKFUN is the intersection between a BALLFUN function and a plane
%   G = DISKFUN(F) is the slice of the BALLFUN F corresponding to 
%   the plane X-Y.
%
%   G = DISKFUN(F, 'x', C) is the slice of the BALLFUN F corresponding to 
%   the plane X = C.
%
%   G = DISKFUN(F, 'y', C) is the slice of the BALLFUN F corresponding to 
%   the plane Y = C.
%
%   G = DISKFUN(F, 'x', C) is the slice of the BALLFUN F corresponding to 
%   the plane Z = C.
%
%   G = DISKFUN(F, PHI, THETA, PSI, C) rotates F using Euler angles phi, theta, 
%   and psi with the ZXZ convention and then evaluates it at the plane Z = C
%   (C = 0 by default)

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
   g = diskfun();
   return
end

% Parse user inputs to get the Euler angles
[phi, theta, psi, c] = parseInputs(varargin{:});

% Rotate f using Euler angles phi, theta and psi
f = rotate(f, phi, theta, psi);

% Get the size
[m,n,~] = size(f);

% If n is odd, make it even
m = m + 1-mod(m,2);
n = n + mod(n,2);

% Evaluation points in [c, 1]
rho = chebpts(m);
rho = rho(ceil(m/2):end);

% Evaluation points in [-pi,pi[
lambda = pi*trigpts(n);

% Build the grid and evaluate at the plane Z = C
r  = sqrt(rho.^2-c^2);

theta = atan2(rho,c);
G = zeros(length(r),n);
for i = 1:length(r)
   G(i,:) = fevalm(f, r(i), lambda, theta(i));
end

% Return the diskfun
g = diskfun(real(G));
end

function [phi, theta, psi, c] = parseInputs(varargin)
% Parse user inputs to DISKFUN.
c = 0;
if nargin == 0
    phi = 0;
    theta = 0;
    psi = 0;
elseif ischar(varargin{1})
    phi = 0;
    if strcmp(varargin{1},'x')
        theta = -pi/2;
        psi = pi/2;
        % Evaluate at Y-Z plane
    elseif strcmp(varargin{1},'y')
        % Evaluate at X-Z plane
        theta = -pi/2;
        psi = 0;
    elseif strcmp(varargin{1},'z')
        % Evaluate at X-Y plane
        theta = 0;
        psi = 0;
    end
    if nargin >= 2
       c = varargin{2}; 
    end
else
    phi = varargin{1};
    theta = 0;
    psi = 0;
    if nargin >= 2
        theta = varargin{2};
    end
    if nargin >= 3
        psi = varargin{3};
    end
    if nargin >= 4
        c = varargin{4};
    end
end
end