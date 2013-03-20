function pass = sampleTest(op, f)
%SAMPLETEST   Test an evaluation of the input OP against a CHEBTECH approximation.
%   SAMPLETEST(OP, F) evaluates both the function OP and its CHEBTECH
%   representation F at one or more points within [-1,1]. The difference of
%   these values is computed, and if this is sufficiently small (relative to
%   F.VSCALE, F.HSCALE, and F.EPSLEVEL) the test passes and returns TRUE. If the
%   difference is large, it returns FALSE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the interpolation points:
n = length(f);
x = f.chebpts(n);

% Set a tolerance:
tol = max(f.epslevel, 1e3*eps) * n * (f.vscale/f.hscale);

% Choose a point to evaluate at:
if ( n == 1 )
    xeval = 0.61; % Pseudo-random test value
else
    % Test a point where the (finite difference) gradient of values is largest:
    [ignored, index] = max(bsxfun(@rdivide, abs(diff(f.values)), diff(x)));
    xeval = ( x(index + 1) + 1.41*x(index) ) / 2.41;
end

% Evaluate the CHEBTECH:
vFun = feval(f, xeval);

% Evaluate the op:
vOp = feval(op, xeval);

% If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
if ( norm(vOp - vFun, inf) > tol )
    pass = false; % :(
else
    pass = true;  % :)
end

end
