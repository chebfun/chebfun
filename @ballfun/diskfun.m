function g = diskfun(f, varargin)
% DISKFUN is the intersection between a BALLFUN function and a plane
%   DISKFUN(F, PHI, THETA, PSI) rotates F using Euler angles phi, theta, 
%   and psi with the ZXZ convention and then evaluate it at the plane Z = 0
%   DISKFUN(f, 'x') is the slice of the BALLFUN function f corresponding to 
%   the plane X = 0.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse user inputs to get the Euler angles
[phi, theta, psi] = parseInputs(varargin{:});

% Rotate f using Euler angles phi, theta and psi
f = rotate(f, phi, theta, psi);

% Get the size
[m,n,~] = size(f);

% If n is odd, make it even
m = m + 1-mod(m,2);
n = n + mod(n,2);

% Evaluation points in [0,1]
r = chebpts(m);
r = r(ceil(m/2):end);

% Evaluation points in [-pi,pi[
lambda = pi*trigpts(n);

% Build the grid and evaluate at the plane Z = 0
[rr, ll, tt] = ndgrid(r, lambda, pi/2);
G = feval(f,rr,ll,tt);

% Return the diskfun
g = diskfun(real(G));
end

function [phi, theta, psi] = parseInputs(varargin)
% Parse user inputs to DISKFUN.
    if nargin == 0
        phi = 0;
        theta = 0;
        psi = 0;
    elseif nargin == 1
        if ischar(varargin{1})
            if strcmp(varargin{1},'x')
                phi = 0;
                theta = -pi/2;
                psi = pi;
                % Evaluate at Y-Z plane
            elseif strcmp(varargin{1},'y')
                % Evaluate at X-Z plane
                phi = 0;
                theta = -pi/2;
                psi = 0;
            elseif strcmp(varargin{1},'z')
                % Evaluate at X-Y plane
                phi = 0;
                theta = 0;
                psi = 0;
            end
        else
            phi = varargin{1};
            theta = 0;
            psi = 0;
        end
    elseif nargin == 2
        phi = varargin{1};
        theta = varargin{2};
        psi = 0;
    else
        phi = varargin{1};
        theta = varargin{2};
        psi = varargin{3};
    end
end