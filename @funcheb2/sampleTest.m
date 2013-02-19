function pass = sampleTest(op, values, vscale, epslevel, pref)
%SAMPLETEST  Test an evaluation of the input OP against the FUNCHEB2 approx.
%   SAMPLETEST(OP, VALUES, VSCALE, EPSLEVEL) evaluates the both the function OP
%   and the interpolant defined by the values VALUES at one or more points
%   within [-1,1] and compares the difference. If this difference is
%   sufficiently small (relative to VSCALE and EPSLEVEL) the the test passes and
%   returns true. If not, it returns false.
%
%   SAMPLETEST(OP, VALUES, VSCALE, EPSLEVEL, PREF) overrides some of the
%   default preferences (in particular HSCALE) with those in the preferences
%   structure PREF.
%
% See also funcheb2.pref.

% Parse the inputs:
if ( nargin < 3 )
    vscale = norm(values(:), inf);
end
if ( nargin < 5)
    pref = funcheb2.pref;
end
if ( nargin < 4 )
    epslevel = pref.funcheb2.eps;
end

% Points of 2nd kind:
n = size(values, 1);
x = funcheb2.chebpts(n); 

% Set a tolerence:
hscale = pref.funcheb2.hscale;
tol = max(epslevel, 1e3*eps) * n * (vscale/hscale);

% Choose a point to evaluate at:
if ( n == 1 )
    xeval = 0.61; % Pseudo-random test value
else
    % Test a point where the (finite difference) gradient of values is largest:
    [ignored, indx] = max( abs(diff(values)) ./ ...
                            repmat(diff(x), 1, size(values, 2)) );
                        % [TODO]: bsxfun?
    xeval = ( x(indx+1) + 1.41*x(indx) ) / 2.41;
end

% Evaluate the FUNCHEB2:
vFun = funcheb2.bary(xeval, values);

% Evaluate the op:
vOp = feval(op, xeval);

% % Adjust size if needed:
% if ( size(vOp, 2) ~= size(vFun, 2) )
%     vOp = repmat(vOp, 1, size(vFun, 2));
% end

% If the FUNCHEB2 evaluation differs from the op evaluation, sampletest failed:
if ( norm( vOp - vFun, inf ) > tol )
    pass = false; % :(
else
    pass = true;  % :)
end

end