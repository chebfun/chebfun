function g = inv2(f,varargin)
%INV2  Invert a CHEBFUN.
%   INV2(F) will attempt to invert the monotonic chebfun F. If F has zero
%   derivatives at its endpoints, then it is advisable to turn Splitting ON.
%
%   INV2(F,'SPLITTING','ON') turns Splitting ON locally for the inv command.
%
%   INV2(F,'EPS',TOL) will construct with the relative tolerance set by TOL. If
%   no tolerance is passed, TOL = 100*chebfunpref('eps') is used. EPS should be
%   set to at least a factor of 100 larger than the accuracy of F.
%
%   INV2(F,'MONOCHECK','ON'/'OFF') turns the check for monotonicity ON or OFF
%   respectively. It is OFF by default.
%
%   G = INV2(F,'RANGECHECK','ON'/'OFF') enforces that the range of G exactly
%   matches the domain of F (by adding a linear function) or not. RANGECHECK OFF
%   is the default behaviour.
%
%   Any of the preferences above can be used in tandem.
%
% Example:
%   f = chebfun(@(x) tanh(7*x)./tanh(7) + 1, [-.5, .5]);
%   g = inv2(f, 'splitting', 'off', 'rangecheck', 'off', 'monocheck', 'off');
%
% Note, this function is experimental and slow! INV may be the better choice for
% piecewise functions, where as INV2 is good for smooth functions.
%
% See also INV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: This method requires a test.

% No quasimatrix support!
if ( min(size(f)) > 1)
    error('CHEBFUN:inv:noquasi','No support for array-valued CHEBFUN objects.');
end

% Default options:
splitYesNo = chebfunpref('splitting');
tol = 100*chebfunpref('eps');
monoCheck = false;
rangeCheck = false;

% Turn splitting on if F is piecewise:
if ( length(f.domain) > 2 && ~splitYesNo )
    splitYesNo = true;
end

% Parse input
while ( numel(varargin) > 1 )
    if ( strcmpi(varargin{1}, 'splitting') )
        splitYesNo = checkOnOff(varargin{2});
    elseif ( strcmpi(varargin{1}, 'eps') )
        tol = varargin{2};
    elseif ( strcmpi(varargin{1}, 'monocheck') )
        monoCheck = checkOnOff(varargin{2});
    elseif ( strcmpi(varargin{1}, 'rangecheck') )
        rangeCheck = checkOnOff(varargin{2});
    else
        error('CHEBFUN:inv2:inputs', ...
            [varargin{1}, 'is an unrecognised input to INV().']);
    end
    varargin(1:2) = [];
end

% Compute the derivative:
fp = diff(f);

% Monotonicity check:
if ( monoCheck )
    tPoints = roots(fp);
    if ( ~isempty(tPoints) )
        endtest = zeros(length(tPoints), 1);
        for k = 1:length(tPoints)
            endtest(k) = min(abs(tPoints(k) - f.domain));
        end
        if ( any(endtest > 100*abs(feval(f, tPoints))*tol) )
            error('CHEBFUN:inv2:notmonotonic', 'chebfun F must be monotonic its domain.');
        elseif ( ~splitYesNo )
            warning('CHEBFUN:inv:singularendpoints', ['F is monotonic, but ', ...
                'INV2(F) has singular endpoints. Suggest you try ''splitting on''.']);
        end
    end
end

% Assign preferences:
p.techPrefs.resampling = 0;
p.enableBreakpointDetection = splitYesNo;
p.techPrefs.eps = tol;
p.techPrefs.minsamples = length(f);
p.techPrefs.sampleTest = 0;

% Compute the inverse:
gDomain = minandmax(f).';
x = chebfun('x', gDomain);
g = chebfun(@(x) fInverse(f, fp, x, tol), gDomain, p );

% Scale so that the range of g is the domain of f:
if ( rangeCheck )
    [gRange, gx] = minandmax(g);
    g = g + (gx(2) - x)*(f.domain(1) - gRange(1))/diff(gx) ...
          + (x - gx(1))*(f.domain(end) - gRange(2))/diff(gx);
end

end

function y = fInverse(f, fp, x, tol)
    tol = tol/5;
    y = zeros(length(x), 1);
    % Newton iteration:
    t = f.domain(1);
    for j = 1:length(x);
        ft = feval(f, t) - x(j);
        while ( abs(ft) > tol ) % Newton step
            fpt = feval(fp, t);
            t = t - ft./fpt;
            ft = feval(f, t) - x(j);
        end
        y(j,1) = t;
    end
end    
    
function value = checkOnOff(value)
    if ( ischar(value) )
        % If ON or OFF used -> change to true or false.
        if ( strcmpi(value,'on') )
            value = true;
        elseif ( strcmpi(value,'off') )
            value = false; 
        end
    else
        value = logical(value);
    end
end
