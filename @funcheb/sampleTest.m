function pass = sampleTest(op, f)
%SAMPLETEST  Test an evaluation of the input OP against a FUNCHEB approximation.
%   SAMPLETEST(OP, F) evaluates the both the function OP and its FUNCHEB
%   representation F at one or more points within [-1,1]. The difference of
%   these values is compared, and if this is is sufficiently small (relative to
%   F.VSCALE, F.HSCALE, and F.EPSLEVEL) the test passes and returns true. If the
%   difference is large, it returns false.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the interpolation points:
n = length(f);
x = f.chebpts(n);

% Set a tolerence:
tol = max(f.epslevel, 1e3*eps) * n * (f.vscale/f.hscale);

% Choose a point to evaluate at:
if ( n == 1 )
    xeval = 0.61; % Pseudo-random test value
else
    % Test a point where the (finite difference) gradient of values is largest:
    [ignored, indx] = max( bsxfun(@rdivide, abs(diff(f.values)), diff(x)) );
    xeval = ( x(indx+1) + 1.41*x(indx) ) / 2.41;
end

% Evaluate the FUNCHEB:
vFun = feval(f, xeval);

% Evaluate the op:
vOp = feval(op, xeval);

% If the FUNCHEB evaluation differs from the op evaluation, sampletest failed:
if ( norm( vOp - vFun, inf ) > tol )
    pass = false; % :(
else
    pass = true;  % :)
end

end
