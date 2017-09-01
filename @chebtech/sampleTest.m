function pass = sampleTest(op, values, f, data, pref)
%SAMPLETEST   Test an evaluation of input OP against a CHEBTECH approximation.
%   SAMPLETEST(OP, VALUES, F) evaluates both the function OP and its CHEBTECH
%   representation F at one or more points within [-1,1]. The difference of
%   these values is computed, and if this is sufficiently small (relative to
%   DATA.VSCALE and DATA.HSCALE) the test passes and returns TRUE. If the
%   difference is large, it returns FALSE. SAMPLETEST(OP, VALUES, F, DATA) will
%   test relative to the values given in DATA.VSCALE, rather than VSCALE(F).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set a tolerance:
tol = sqrt(max(eps, pref.chebfuneps));

% Scale TOL by the MAX(DATA.HSCALE*||F||, DATA.VSCALE). 
% (See standardCheck for explanation)
vscaleF = max(abs(values), [], 1);
tol = tol.*max(data.hscale*vscaleF, data.vscale);

% choose points to evaluate
xeval = [-0.357998918959666; 0.036785641195074];

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
