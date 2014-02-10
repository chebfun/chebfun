function out = isdecay(f)
%ISDECAY   Test if a CHEBTECH decays faster than a single root at endpoints.
%   ISDECAY(F) returns a 1x2 row vector each of which indicates whether F
%   vanishes at one of the endpoints faster than a sinhle root. An entry TRUE is
%   returned if F has a boundary root with multiplicity larger than one, FALSE
%   otherwise. 
%
%   Note that ISDECAY is designed and expected to be called only by UNBNDFUN
%   class for handling functions defined on unbounded domains.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

out = zeros(1, 2);

% Set a tolerance:
tol = 1e3*get(f, 'epslevel')*get(f, 'vscale');

% Get the points at which the function is sampled and the function values:
x = get(f, 'points');

%% Right endpoint:

% Copy F:
g = f;

% TODO: This should be done in coefficient space (more stable).
% Peel off a factor of (1-x):
g.values = g.values./(1-x);

% Extrapolate the boundary value:
g.values = extrapolate(g);

% Decaying or growing?
rate = diff(abs(g.values(end-1:end)));

% Check decaying speed:
if all( abs(g.values(end-4:end)) < 1e2*tol ) || ... % If exponentially decay
   ( ( abs(g.values(end)) < tol ) && ( rate < 0 ) ) % If decays fast enough
    out(2) = 1;
end

%% Left endpoint:

% Copy F:
g = f;

% TODO: This should be done in coefficient space (more stable).
% Peel off a factor of (1+x):
g.values = g.values./(1+x);

% Extrapolate the boundary value:
g.values = extrapolate(g);

% Decaying or growing?
rate = diff(abs(g.values(1:2)));

% Check decaying speed:
if all( abs(g.values(1:5)) < 1e2*tol ) || ... % If exponentially decay
   ( ( abs(g.values(1)) < tol ) && ( rate > 0 ) ) % If decays fast enough
    out(1) = 1;
end

end