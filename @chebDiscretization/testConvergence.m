function [isDone, epsLevel, vscale] = testConvergence(disc, values, vscale, pref)
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
%      EPSLEVEL: Apparent resolution accuracy (relative to VSCALE or the
%      functions' intrinsic scale).

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Document inputs.

if ( nargin < 4 )
    pref = cheboppref;
    if ( nargin < 3 )
        vscale = 0;   % will have no effect
    end
end

% We will test on an arbitrary linear combination of the individual functions.
s = 1 ./ (3*(1:numel(values))).';
newValues = cell2mat(values(:).')*s;

% Convert to a piecewise CHEBFUN.
u = toFunctionOut(disc, newValues);

% Test convergence on each piece. Start by obtaining the Chebyshev coefficients
% of all pieces, which we can then pass down to the testPiece method:
coeffs = get(u, 'coeffs');
d = disc.domain;
numInt = numel(d) - 1;
isDone = false(1, numInt);
epsLevel = 0;

% If an external vscale was supplied, it can supplant the inherent scale of the
% result.
vscale = max(u.vscale, vscale);
prefTech = chebtech.techPref();
prefTech.eps = pref.errTol;

for i = 1:numInt
    f = chebtech2( {[],coeffs{i}} );
    f.vscale = vscale;
    [isDone(i), neweps] = plateauCheck(f, newValues, prefTech);
    epsLevel = max(epsLevel, neweps);
end

end