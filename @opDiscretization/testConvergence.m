function [isDone, cutoff, vscale] = testConvergence(disc, values, vscale, pref)
%TESTCONVERGENCE   Happiness check.
%   Given: 
%      DISC: chebDiscretization, 
%      VALUES: a cell array of scalars/sampled function values (see the
%           toFunction method),
%      VSCALE: scalar value giving the desired scale by which to measure
%           relative convergence against (defaults to 0, which means use
%           the intrinsic scale of the result),
%      PREF: A cheboppref() options structure.
%
%   Output:  
%      ISDONE: True if the functions passed in are sufficiently resolved.
%      CUTOFF: The point at which the coefficient series for DISC should
%              be chopped. If ~ISDONE then CUTOFF = LENGTH(DISC).
%      VSCALE: Maximum of the input VSCALE and the computed VSCALE of DISC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    vscale = 0;   % will have no effect
end
if ( nargin < 4 )
    pref = cheboppref;
end

% Convert to a piecewise array-valued CHEBFUN.
u = toFunctionOut(disc, cat(2, values{:}));
numCol = size(u, 2);

% This is a cell array of coefficients (one for each piece).
coeffs = get(u, 'coeffs', 1);

d = disc.domain;
numInt = numel(d) - 1;
isDone = false(numInt, 1);
cutoff = zeros(numInt, numCol);

% Get the discretization, and the appropriate tech to use:
discPreference = pref.discretization();
tech = discPreference.returnTech();
tech = tech();

% If an external vscale was supplied, it can supplant the inherent scale of the
% result.
vscale = max(u.vscale, vscale);
data.vscale = vscale;
% TODO:  Assign hscale (to data.hscale)?
prefTech = tech.techPref();
prefTech.chebfuneps = pref.bvpTol;
prefTech.happinessCheck = pref.happinessCheck;

for i = 1:numInt
    c = cat(2, coeffs{i,:});
    f = tech.make({[], c});
    [isDone(i), cutoff(i,:)] = happinessCheck(f, [], [], data, prefTech);
end

isDone = all(isDone, 2);

end
