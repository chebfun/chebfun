function g = inv(f, varargin)
%INV   Invert a CHEBFUN.
%   INV(F) will attempt to invert the monotonic chebfun F. If F has zero
%   derivatives at its endpoints, then it is advisable to turn Splitting ON.
%
%   INV(F, 'SPLITTING', 'ON') turns Splitting ON locally for the inv command.
%
%   INV(F, 'EPS', TOL) will construct with the relative tolerance set by TOL. If
%   no tolerance is passed, TOL = chebfunpref('eps') is used.
%
%   INV(F, 'MONOCHECK', 'ON'/'OFF') turns the check for monotonicity ON or OFF
%   respectively. It is OFF by default.
%
%   G = INV(F, 'RANGECHECK', 'ON'/'OFF') enforces that the range of G exactly
%   matches the domain of F (by adding a linear function) or not. RANGECHECK OFF
%   is the default behaviour.
%
%   Any of the preferences above can be used in tandem.
%
% Example: 
%   x = chebfun('x');
%   f = sign(x) + x;
%   g = inv(f, 'splitting', true);
%
%   Note, this function is experimental and slow! INV may be the better
%   choice for piecewise functions, where as INV2 is good for smooth
% f unctions.
%
% See also INV2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% No quasimatrix support:
if ( min(size(f)) > 1 )
    error('CHEBFUN:inv:noquasi', 'No support for array0valued CHEBFUN objects.');
end

% Default options:
p = chebpref;
splitYesNo = p.enableSingularityDetection;
tol = epslevel(f);
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
        error('CHEBFUN:inv:inputs', ...
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
            error('CHEBFUN:inv:notmonotonic', 'chebfun F must be monotonic its domain.');
        elseif ( ~splitYesNo )
            warning('CHEBFUN:inv:singularendpoints', ['F is monotonic, but ', ...
                'INV(F) has singular endpoints. Suggest you try ''splitting on''.']);
        end
    end
end

% Compute the inverse:
gDomain = minandmax(f).';
x = chebfun(@(x) x, gDomain);
g = chebfun(@(x) fInverse(f, x), gDomain, 'resampling', 0, 'splitting', splitYesNo, 'eps', tol);

% Scale so that the range of g is the domain of f:
if ( rangeCheck )
    [gRange, gx] = minandmax(g);
    g = g + (gx(2)-x)*(f.domain(1) - gRange(1))/diff(gx) ...
          + (x-gx(1))*(f.domain(end) - gRange(2))/diff(gx);
end
  
    function y = fInverse(f, x)
        % FINVERSE(F, X) attemptes to return the values of F^{-1}(X) using
        % CHEBFUN/ROOTS().
        y = zeros(length(x), 1);
        % Vectorise:
        for j = 1:length(x)
            temp = roots(f - x(j));
            if ( length(temp) ~= 1 )
                fvals = feval(f, f.domain);
                err = abs(fvals - x(j));
                [temp, k] = min(err);
                if ( err(k) > 100*tol*abs(fvals(k)))
                    error('CHEBFUN:inv:notmonotonic2', 'Chebfun must be monotonic.');
                end
            end
            y(j, 1) = temp;
        end
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


