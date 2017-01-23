function g = inv(f, varargin)
%INV   Invert a CHEBFUN.
%   FINV = INV(F) attempts to compute the inverse of the monotonic CHEBFUN F.
%
%   FINV = INV(..., 'ALGORITHM', ALGSTR) selects the algorithm used to compute
%   the values of the inverse of F.  Possible values for ALGSTR are:
%      'ROOTS'  - Compute the inverse using ROOTS().
%      'NEWTON' - Compute the inverse using a Newton iteration.
%      'BISECTION' - Compute the inverse using bisection as the rootfinder.
%      'REGULAFALSI' - Compute the inverse using Regula Falsi as the rootfinder.
%      'ILLINOIS' - Compute the inverse using Illinois as the rootfinder.
%      'BRENT' - Compute the inverse using Brent's method as the rootfinder.
%   The default algorithm is 'BRENT'.
%
%   FINV = INV(F, PREF) uses the preferences specified by the structure or
%   CHEBFUNPREF object PREF when constructing the inverse.
%
%   FINV = INV(..., 'EPS', TOL) will construct with the relative tolerance set
%   by TOL.  If no tolerance is passed, TOL = EPS is used.
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
%      f = x + .5*abs(x) + .6*sign(x-.5);
%      g = inv(f);
%      plot(f, x, 'b', g, '--r') % <-- Note, plot(f, x) not plot(x, f).
%
% See also ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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

% Find the domain of the output:
gEnds = minandmax(f).';
fBreaks = f.domain(2:end-1);
gBreaksL = feval(f, fBreaks, 'left');  % We must evaluate to the left and the 
gBreaksR = feval(f, fBreaks, 'right'); % of breaks in f in case there are jumps.
gBreaks = chebfun.tolUnion(gBreaksL, gBreaksR);
gDomain = chebfun.tolUnion(gEnds, gBreaks);

% Compute the inverse:
if ( opts.algorithm == 1 )     % Algorithm based on ROOTS.
    g = chebfun(@(x) fInverseRoots(f, x, tol), gDomain, pref, 'noVectorCheck');
elseif ( opts.algorithm == 2 ) % Newton iteration algorithm.
    g = chebfun(@(x) fInverseNewton(f, fp, x, tol), gDomain, pref, 'noVectorCheck');
elseif ( opts.algorithm == 3 ) % Bisection based algorithm.
    g = chebfun(@(x) fInverseBisection(f, x), gDomain, pref, 'noVectorCheck');
elseif ( opts.algorithm == 4 ) % Regula Falsi based algorithm.
    g = chebfun(@(x) fInverseRegulaFalsi(f, x), gDomain, pref, 'noVectorCheck');
elseif ( opts.algorithm == 5 ) % Illinois based algorithm.
    g = chebfun(@(x) fInverseIllinois(f, x), gDomain, pref, 'noVectorCheck');
elseif ( opts.algorithm == 6 ) % Brent's method.
    g = chebfun(@(x) fInverseBrent(f, x), gDomain, pref, 'noVectorCheck');    
else
    error('CHEBFUN:CHEBFUN:inv:algorithm', 'Invalid algorithm selected.');
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
tol = eps;
opts.monoCheck = false;
opts.rangeCheck = false;
opts.algorithm = 6; % Default  = brent's method.

% Parse preference input:
if ( (nargin > 1) && isa(varargin{1}, 'chebfunpref') )
    pref = varargin{1};
    varargin(1) = [];
else
    pref = chebfunpref();
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
        elseif ( strcmpi(varargin{2}, 'regulafalsi') )
            opts.algorithm = 4;            
        elseif ( strcmpi(varargin{2}, 'illinois') )
            opts.algorithm = 5;
        elseif ( strcmpi(varargin{2}, 'brent') )
            opts.algorithm = 6;            
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
if ( opts.algorithm == 2 );
    % TODO:  CHEBFUN is not supposed to set the refinementFunction preference
    % because it doesn't belong to the list of "abstract" preferences required
    % of all techs.  Do we really need to alter it here?
    pref.techPrefs.refinementFunction = 'resampling';
end
pref.techPrefs.chebfuneps = tol;
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

a = f.domain(1);
b = f.domain(end);
c = (a + b)/2;

% The loop below is written for functions which are monotone increasing.
% Flip the signs if this is not the case.
sgn = sign(diff(feval(f, [a b])));

while ( norm(b - a, inf) >= eps )   
    vals = sgn*(feval(f, c) - x);
    % Bisection:
    I1 = (vals <= -eps);
    I2 = (vals >=  eps);
    I3 = ~I1 & ~I2;
    a = I1.*c + I2.*a + I3.*c;
    b = I1.*b + I2.*c + I3.*c;
    c = (a + b)/2;
end
y = c;

end

function y = fInverseRegulaFalsi(f, x)
%FINVERSEREGULAFALSI(F, X)   Compute F^{-1}(X) using Regula Falsi.
a = f.domain(1);
b = f.domain(end);
fa = feval(f, a) - x;
fb = feval(f, b) - x;
c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cOld = inf;

while ( norm(c - cOld, inf) >= eps )   
    cOld = c;
    fc = feval(f, c) - x;
    
    I1 = (fc < 0);
    I2 = (fc > 0);
    I3 = ~I1 & ~I2;
    a = I1.*c + I2.*a + I3.*c;
    b = I1.*b + I2.*c + I3.*c;
    fa = I1.*fc + I2.*fa + I3.*fc;
    fb = I1.*fb + I2.*fc + I3.*fc;
    step = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
    
end
y = c;

end

function y = fInverseIllinois(f, x)
%FINVERSEILLINOIS(F, X)   Compute F^{-1}(X) using the Illinois algorithm.
a = f.domain(1);
b = f.domain(end);
fa = feval(f, a) - x;
fb = feval(f, b) - x;
c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cOld = inf;

side = zeros(size(x));

while ( norm(c - cOld, inf) >= eps )   
    cOld = c;
    fc = feval(f, c) - x;
    
    I1 = (fc < 0);
    I2 = (fc > 0);
    I3 = ~I1 & ~I2;
    a = I1.*c + I2.*a + I3.*c;
    b = I1.*b + I2.*c + I3.*c;
    fa = I1.*fc + I2.*fa + I3.*fc;
    fb = I1.*fb + I2.*fc + I3.*fc;
    
    fb(side(I1) == -1) = fb(side(I1) == -1)/2;
    side(I1) = -1;
    fa(side(I1) == 1) = fa(side(I1) == 1)/2;
    side(I2) = 1;
    
    step = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
    
end
y = c;

end

function y = fInverseBrent(f, x)
%FINVERSEBRENT(F, X)   Compute F^{-1}(X) using Brent's method.

% References:
% [1] en.wikipedia.org/w/index.php?title=Brent%27s_method&oldid=610174347
% [2] Brent, R. P. (1973), "Chapter 4: An Algorithm with Guaranteed Convergence 
%     for Finding a Zero of a Function", Algorithms for Minimization without 
%     Derivatives, Englewood Cliffs, NJ: Prentice-Hall, (1973).

% Set a and b:
a = f.domain(1);
b = f.domain(end);

% Calculate f(a) and f(b) (including x shift):
fa = feval(f, a) - x;
fb = feval(f, b) - x;
fs = inf;

% Make a and b vectors:
z = zeros(size(x));
a = a + z;
b = b + z;

% if |f(a)| < |f(b)| then swap (a,b) end if:
idx = abs(fa) < abs(fb);

tmp = a(idx);
a(idx) = b(idx);
b(idx) = tmp;

tmp = fa(idx);
fa(idx) = fb(idx);
fb(idx) = tmp;

% Set c = a;
c = a;
fc = fa;

% Set mFlag and delta:
mFlag = true(size(x)); % <-- Stores previous decision state.
delt = eps;

% Initialise these too:
s = a;
d = c;

while ( (norm(fs, inf) > eps) && ...
        (norm(b - a, inf) > eps*norm([b(:) ; a(:)], inf))  )
    % Inverse quadratic:
    s_iq = a.*fb.*fc./((fa-fb).*(fa-fc)) + b.*fa.*fc./((fb-fa).*(fb-fc)) + ...
        c.*fa.*fb./((fc-fa).*(fc-fb));
    % Secant:
    s_sc = b - fb.*(b-a)./(fb-fa);
    % Bisection:
    s_bi = (a + b)/2;
    
    % Decide which update to use. Essentially the same as described on the
    % Wikipedia page given above (permalink).
    
    % if f(a) ~= f(c) and f(b) ~= f(c) then (vectorized)
    idx = (fa ~= fc) & (fb ~= fc);
    s(idx)  = s_iq(idx);  % Take the inverse quadratic step.
    s(~idx) = s_sc(~idx); % Take the secant step.

    % Conditions for bisection:
    idx = ( (3*a+b)/4 < b  & (s < (3*a+b)/4 | s > b) ) | ...     % condition 1a
          ( b <= (3*a+b)/4 & (s < b | s > (3*a+b)/4) ) | ...     % condition 1b
          ( mFlag  & abs(s-b) >= abs(b-c)/2 ) | ...              % condition 2
          ( ~mFlag & abs(s-b) >= abs(c-d)/2 ) | ...              % condition 3
          ( mFlag  & abs(b-c) < delt ) | ...                     % condition 4
          ( ~mFlag & abs(c-d) < delt ) | ...                     % condition 5
          ( abs(b-a) < delt );
    s(idx) = s_bi(idx);   % Take the bisection step.
    mFlag = idx;  % <-- Stores previous decision state.
    
    % Calculate f(s) (including x shift):
    fs = feval(f, s) - x;
    
    % Store as previousol values for next iteration"
    d = c;
    c = b;
    fc = fb;
    
    % if f(a) f(s) < 0 then b := s else a := s end if (vectorized)
    idx = fa.*fs <= 0;
    b(idx) = s(idx);
    fb(idx) = fs(idx);
    a(~idx) = s(~idx);
    fa(~idx) = fs(~idx);

    % if |f(a)| < |f(b)| then swap (a,b) end if: (vectorized)
    idx = abs(fa) < abs(fb);
    tmp = a(idx);
    a(idx) = b(idx);
    b(idx) = tmp;
    tmp = fa(idx);
    fa(idx) = fb(idx);
    fb(idx) = tmp;

end

% Output y:
y = s;

end
