function g = inv(f, varargin)
%INV   Invert a CHEBFUN.
%   FINV = INV(F) attempts to compute the inverse of the monotonic CHEBFUN F.
%
%   FINV = INV(..., 'ALGORITHM', ALGSTR) selects the algorithm used to compute
%   the values of the inverse of F.  Possible values for ALGSTR are:
%      'ROOTS'  - Compute the inverse using ROOTS().
%      'NEWTON' - Compute the inverse using a Newton iteration.
%      'BISECTION' - Compute the inverse using bisection as the rootfinder.
%   The default algorithm is 'BISECTION'.
%
%   FINV = INV(F, PREF) uses the preferences specified by the structure or
%   CHEBFUNPREF object PREF when constructing the inverse.
%
%   FINV = INV(..., 'EPS', TOL) will construct with the relative tolerance set
%   by TOL.  If no tolerance is passed, TOL = EPSLEVEL(F) is used.
%
%   FINV = INV(..., 'SPLITTING', 'ON') enables breakpoint detection locally for
%   INV.  Setting this option (or the equivalent preference in PREF) is
%   particularly advisable if F has zero derivatives at its endpoints.
%
%   FINV = INV(..., 'MONOCHECK', 'ON') enables an optional check for
%   monotonicity.
%
%   FINV = INV(..., 'RANGECHECK', 'ON') enforces that the range of FINV exactly
%   matches the domain of F (by adding a linear function).
%
%   Any of the name-value option pairs listed above can be used in tandem.
%
%   Example:
%      x = chebfun('x');
%      f = sign(x) + x;
%      g = inv(f, 'splitting', 'on');
%
%   NB:  This function is experimental and slow!  Use of the 'BISECTION'
%   (default) and 'ROOTS' algorithm may be the better choice for piecewise
%   functions, whereas the 'NEWTON' algorithm is good for smooth functions.
%
% See also ROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% No quasimatrix support:
if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:inv:noquasi', ...
        'INV does not support array-valued CHEBFUN objects or quasimatrices.');
end

% Parse the inputs:
[tol, opts, pref] = parseInputs(f, varargin{:});

if ( opts.monoCheck )
    % Compute the derivative:
    fp = diff(f);
    % Monotonicity check:
    doMonoCheck(f, fp, tol, pref.splitting);
else
    fp = [];
end

% Compute the inverse:
gDomain = minandmax(f).';
if ( opts.algorithm == 1 )     % Algorithm based on ROOTS.
    g = chebfun(@(x) fInverseRoots(f, x, tol), gDomain, pref);
elseif ( opts.algorithm == 2 ) % Newton iteration algorithm.
    g = chebfun(@(x) fInverseNewton(f, fp, x, tol), gDomain, pref);
elseif ( opts.algorithm == 3 ) % Bisection based algorithm.
    g = chebfun(@(x) fInverseBisection(f, x), gDomain, pref);
end

% Scale so that the range of g is the domain of f:
if ( opts.rangeCheck )
    g = adjustRange(g, f.domain);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [tol, opts, pref] = parseInputs(f, varargin)

% Default options:
tol = epslevel(f);
opts.monoCheck = false;
opts.rangeCheck = false;
opts.algorithm = 3;

% Parse preference input:
if ( (nargin > 1) && isa(varargin{1}, 'chebfunpref') )
    pref = varargin{1};
    varargin(1) = [];
else
    pref = chebfunpref();
end

% Enable breakpoint detection if F is piecewise:
if ( length(f.domain) > 2 )
    p.splitting = true;
end

% Parse name/value pairs.
while ( numel(varargin) > 1 )
    if ( strcmpi(varargin{1}, 'splitting') )
        pref.splitting = checkOnOff(varargin{2});
    elseif ( strcmpi(varargin{1}, 'eps') )
        tol = varargin{2};
    elseif ( strcmpi(varargin{1}, 'monocheck') )
        opts.monoCheck = checkOnOff(varargin{2});
    elseif ( strcmpi(varargin{1}, 'rangecheck') )
        opts.rangeCheck = checkOnOff(varargin{2});
    elseif ( strcmpi(varargin{1}, 'algorithm') )
        if ( strcmpi(varargin{2}, 'roots') )
            opts.algorithm = 1;
        elseif ( strcmpi(varargin{2}, 'newton') )
            opts.algorithm = 2;
        elseif ( strcmpi(varargin{2}, 'bisection') )
            opts.algorithm = 3;
        else
            error('CHEBFUN:CHEBFUN:inv:badAlgo', ...
                'Unrecognized value for ''algorithm'' input.');
        end
    else
        error('CHEBFUN:CHEBFUN:inv:inputs', ...
            [varargin{1} ' is an unrecognised input to INV().']);
    end
    varargin(1:2) = [];
end

% Assign preferences:
pref.techPrefs.resampling = 1;
pref.techPrefs.eps = tol;
pref.techPrefs.minSamples = length(f);
pref.techPrefs.sampleTest = 0;

end

function value = checkOnOff(value)
%CHECKONOFF   Convert 'on' and 'off' strings to their logical equivalents.

if ( ischar(value) )
    if ( strcmpi(value,'on') )
        value = true;
    elseif ( strcmpi(value,'off') )
        value = false;
    end
else
    value = logical(value);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monotonicity and range checks.

function doMonoCheck(f, fp, tol, splitYesNo)

tPoints = roots(fp);
if ( ~isempty(tPoints) )
    endtest = zeros(length(tPoints), 1);
    for k = 1:length(tPoints)
        endtest(k) = min(abs(tPoints(k) - f.domain));
    end
    if ( any(endtest > 100*abs(feval(f, tPoints))*tol) )
        error('CHEBFUN:CHEBFUN:inv:doMonoCheck:notMonotonic', ...
            'F must be monotonic on its domain.');
    elseif ( ~splitYesNo )
        warning('CHEBFUN:CHEBFUN:inv:doMonoCheck:singularEndpoints', ...
            ['F is monotonic, but its inverse has singular endpoints. ' ...
            'Enabling breakpoint detection is advised.']);
    end
end

end

function g = adjustRange(g, fDomain)

x = chebfun(@(x) x, g.domain);
[gRange, gx] = minandmax(g);
g = g + (gx(2) - x)*(fDomain(1) - gRange(1))/diff(gx) ...
    + (x - gx(1))*(fDomain(end) - gRange(2))/diff(gx);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse finding algorithms.

function y = fInverseRoots(f, x, tol)
%FINVERSEROOTS(F, X, TOL)   Compute F^{-1}(X) using CHEBFUN.ROOTS().

y = zeros(length(x), 1);
% Vectorise:
for j = 1:length(x)
    temp = roots(f - x(j));
    if ( length(temp) ~= 1 )
        fvals = feval(f, f.domain);
        err = abs(fvals - x(j));
        [ignored, k] = min(err);
        if ( err(k) > 100*tol*abs(fvals(k)))
            error('CHEBFUN:CHEBFUN:inv:notmonotonic2', 'f must be monotonic.');
        end
        temp = temp(k);
    end
    y(j, 1) = temp;
end

end

function y = fInverseNewton(f, fp, x, tol)
%FINVERSENEWTON(F, X)   Compute F^{-1}(X) using a Newton iteration.

tol = tol/5;
y = zeros(length(x), 1);

% compute the derivative of f:
if ( isempty(fp) )
    fp = diff(f);
end

% If the sample points are dense enough in the construction domain that our
% initial guesses scheme for the Newton iteration will work (see below), use
% that.  Otherwise, compute the inverse with ROOTS.
if ( length(x) >= length(f) )
    t = f.domain(1);
    for j = 1:length(x);
        % Do Newton iteration with solution for previous sample point as the
        % starting guess:
        ft = feval(f, t) - x(j);
        counter = 0;
        while ( abs(ft) > tol )
            fpt = feval(fp, t);
            t = t - ft./fpt;
            ft = feval(f, t) - x(j);
            counter = counter + 1;
            if ( counter > 10 )
                % TODO:  Issue a warning here?
                break;
            end
        end
        y(j,1) = t;
    end
else
    y = fInverseRoots(f, x, tol);
end

end

function y = fInverseBisection(f, x)
%FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.

a = f.domain(1)*ones(length(x), 1);
b = f.domain(end)*ones(length(x), 1);
c = (a + b)/2;

while ( norm(b - a, inf) >= eps )   
    vals = feval(f, c);
    % Bisection:
    I1 = ((vals-x) <= -eps);
    I2 = ((vals-x) >= eps);
    I3 = ~I1 & ~I2;
    a = I1.*c + I2.*a + I3.*c;
    b = I1.*b + I2.*c + I3.*c;
    c = (a+b)/2;
end

y = c;

end
