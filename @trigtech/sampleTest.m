function pass = sampleTest(op, f, pref)
%SAMPLETEST   Test an evaluation of input OP against a TRIGTECH approximation.
%   SAMPLETEST(OP, F) evaluates both the function OP and its TRIGTECH
%   representation F at one or more points within [-1,1). The difference of
%   these values is computed, and if this is sufficiently small (relative to
%   VSCALE(F)) the test passes and returns TRUE. If the difference is large,
%   it returns FALSE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the interpolation points:
n = length(f);

% Set a tolerance:
tol =sqrt(max(pref.chebfuneps, eps));
tol = tol.*max(vscale(f));

% pseudo random sample points
xeval = [-0.357998918959666; 0.036785641195074];

% Evaluate the TRIGTECH:
vFun = feval(f, xeval);

% Evaluate the op:
vOp = feval(op, xeval);

% If the TRIGTECH evaluation differs from the op evaluation, SAMPLETEST failed:
err = abs(vOp - vFun); % Relative (to vscale) error.
if ( all(max(err) <= tol) )
    pass = true;  % :)
else
    pass = false; % :(
end

end
