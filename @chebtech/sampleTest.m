function pass = sampleTest(op, values, f, vscl)
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
tol = max(f.epslevel, 1e3*eps) * n;

if ( nargin < 4 || isempty(vscl) )
    vscl = max(abs(values), [], 1);
end

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
err = bsxfun(@rdivide, abs(vOp - vFun), vscl); % Relative (to vscl) error.
if ( any(max(abs(err)) > tol) )
    pass = false; % :(
else
    pass = true;  % :)
end

end
