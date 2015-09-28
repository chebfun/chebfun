function pass = sampleTest(op, values, f, vscl, pref)
%SAMPLETEST   Test an evaluation of input OP against a CHEBTECH approximation.
%   SAMPLETEST(OP, VALUES, F) evaluates both the function OP and its CHEBTECH
%   representation F at one or more points within [-1,1]. The difference of
%   these values is computed, and if this is sufficiently small (relative to
%   F.VSCALE, F.HSCALE, and F.EPSLEVEL) the test passes and returns TRUE. If
%   the difference is large, it returns FALSE. SAMPLETEST(OP, VALUES, F, VSCL)
%   will test relative the the values given in VSCL, rather than F.VSCALE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: Describe where we evaluate? (Approx to largest derivative and at
% -1+1e-12, 1-1e-12?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the interpolation points:
n = length(f);
x = f.chebpts(n);

% Set a tolerance:
tol = max(max(f.epslevel, pref.eps), 1e3*pref.sampleTestEps) * n;

if ( nargin < 4 || isempty(vscl) )
    vscl = max(abs(values), [], 1);
end

% Scale TOL by the MAX(F.HSCALE, VSCL/||F||).
% This choice of scaling is the result of undesirable behavior when using
% standardCheck to construct the function f(x) = sqrt(1-x) on the interval [0,1]
% with splitting turned on. Due to the way standardChop checks for plateaus, the
% approximations on the subdomains were chopped incorrectly leading to poor
% quality results. This choice of scaling corrects this by giving less weight to
% subintervals that are much smaller than the global approximation domain, i.e.
% HSCALE >> 1. For functions on a single domain with no breaks, this scaling has
% no effect since HSCALE = 1. 
nrmf = max(abs(values), [], 1);
if ( isempty(vscl) )
    vscl = nrmf;
end
tol = tol.*max(f.hscale*nrmf, vscl);

% Choose a point to evaluate at:
if ( n == 1 )
    xeval = 0.61; % Pseudo-random test value
else
    % Test a point where the (finite difference) gradient of values is largest:
    [ignored, index] = max(bsxfun(@rdivide, abs(diff(values)), diff(x)));
    xeval = ( x(index + 1) + 1.41*x(index) ) / 2.41;
end
xeval = [-1+1e-12 ; xeval ; 1-1e-12];

% Evaluate the CHEBTECH:
vFun = feval(f, xeval);

% Evaluate the op:
vOp = feval(op, xeval);

% If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
err = abs(vOp - vFun); % Relative (to vscl) error.
if ( all(max(abs(err)) <= tol) )
    pass = true;  % :)
else
    pass = false; % :(
end

end
