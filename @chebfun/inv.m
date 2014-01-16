function g = inv(f, varargin)
%INV   Invert a CHEBFUN.
%   FINV = INV(F) attempts to compute the inverse of the monotonic CHEBFUN F.
%
%   FINV =INV(F, PREF) uses the preferences specified by the structure or
%   CHEBPREF object PREF when constructing the inverse.
%
%   FINV = INV(..., 'SPLITTING', 'ON') enables breakpoint detection locally for
%   INV.  Setting this option (or the equivalent preference in PREF) is
%   particularly advisable if F has zero derivatives at its endpoints.
%
%   FINV = INV(..., 'EPS', TOL) will construct with the relative tolerance set
%   by TOL.  If no tolerance is passed, TOL = EPSLEVEL(F) is used.
%
%   FINV = INV(..., 'MONOCHECK', 'ON') enables an optional check for
%   monotonicity.
%
%   FINV = INV(..., 'RANGECHECK', 'ON') enforces that the range of FINV exactly
%   matches the domain of F (by adding a linear function).
%
%   FINV = INV(..., 'ALGORITHM', ALGSTR) selects the algorithm used to compute
%   the values of the inverse of F.  Possible values for ALGSTR are:
%      'ROOTS'  - Compute the inverse using ROOTS().
%      'NEWTON' - Compute the inverse using a Newton iteration.
%   The default algorithm is 'ROOTS'.
%
%   Any of the name-value option pairs listed above can be used in tandem.
%
%   Example:
%      x = chebfun('x');
%      f = sign(x) + x;
%      g = inv(f, 'splitting', 'on');
%
%   NB:  This function is experimental and slow!  Use of the 'ROOTS' algorithm
%   (default) may be the better choice for piecewise functions, whereas the
%   'NEWTON' algorithm is good for smooth functions.
%
% See also ROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% No quasimatrix support:
if ( min(size(f)) > 1 )
    error('CHEBFUN:inv:noquasi', ...
        'INV does not support array-valued CHEBFUN objects or quasimatrices.');
end

% Parse the inputs:
[tol, opts, pref] = parseInputs(f, varargin{:});

% Compute the derivative:
fp = diff(f);

% Monotonicity check:
if ( opts.monoCheck )
    doMonoCheck(f, fp, tol, pref.enableBreakpointDetection);
end

% Compute the inverse:
gDomain = minandmax(f).';
if ( opts.algorithm == 1 )     % Algorithm based on ROOTS.
    g = chebfun(@(x) fInverseRoots(f, x), gDomain, pref);
elseif ( opts.algorithm == 2 ) % Newton iteration algorithm.
    g = chebfun(@(x) fInverseNewton(f, fp, x, tol), gDomain, pref);
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
opts.algorithm = 1;

% Parse preference input:
if ( (nargin > 1) && isa(varargin{1}, 'chebpref') )
    pref = varargin{1};
    varargin(1) = [];
else
    pref = chebpref();
end

% Enable breakpoint detection if F is piecewise:
if ( length(f.domain) > 2 )
    p.enableBreakpointDetection = true;
end

% Parse name/value pairs.
while ( numel(varargin) > 1 )
    if ( strcmpi(varargin{1}, 'splitting') )
        pref.enableBreakpointDetection = checkOnOff(varargin{2});
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
        else
            error('CHEBFUN:inv:badAlgo', ...
                'Unrecognized value for ''algorithm'' input.');
        end
    else
        error('CHEBFUN:inv:inputs', ...
            [varargin{1} ' is an unrecognised input to INV().']);
    end
    varargin(1:2) = [];
end

% Assign preferences:
pref.techPrefs.resampling = 1;
pref.techPrefs.eps = tol;
pref.techPrefs.minsamples = length(f);
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
        error('CHEBFUN:inv:notMonotonic', ...
            'F must be monotonic on its domain.');
    elseif ( ~splitYesNo )
        warning('CHEBFUN:inv:singularEndpoints', ...
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

function y = fInverseRoots(f, x)
%FINVERSEROOTS(F, X)   Compute F^{-1}(X) using CHEBFUN.ROOTS().

y = zeros(length(x), 1);
% Vectorise:
for j = 1:length(x)
    temp = roots(f - x(j));
    if ( length(temp) ~= 1 )
        fvals = feval(f, f.domain);
        err = abs(fvals - x(j));
        [temp, k] = min(err);
        if ( err(k) > 100*tol*abs(fvals(k)))
            error('CHEBFUN:inv:notmonotonic2', 'f must be monotonic.');
        end
    end
    y(j, 1) = temp;
end

end

function y = fInverseNewton(f, fp, x, tol)
%FINVERSENEWTON(F, X)   Compute F^{-1}(X) using a Newton iteration.

tol = tol/5;
y = zeros(length(x), 1);

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
    y = fInverseRoots(f, x);
end

end
