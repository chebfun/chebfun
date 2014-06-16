function d = dirac(f, varargin)
%DIRAC    Dirac delta function.
% D = DIRAC(F) returns a CHEBFUN D which is zero on the domain of the CHEBFUN F
% except at the simple roots of F, where it is infinite.
%
% DIRAC(F, N) is the nth derivative of DIRAC(F).
%
% DIRAC(F) is not defined if F has a zero of order greater than one within the
% domain of F.
%
% If F has break-points, they should not coincide with the roots of F. However,
% F can have simple roots at either end points of its domain.
%
% See also HEAVISIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty argument:
if ( isempty(f) )
    d = f;
    return
end

% Deal with quasimatrices:
if ( numColumns(f) > 1 )
    f = cheb2cell(f);
    for k = 1:numel(f)
        f{k} = dirac(f{k}, varargin{:});
    end
    d = horzcat(f{:});
    return
end

% Handle the case for derivatives of delta function:
if ( nargin > 1 )
    if ( nargin > 2 )
        error('CHEBFUN:CHEBFUN:dirac:dirac', 'Too many input arguments.');
    end
    
    % Order of the derivative of dirac delta function:
    n = varargin{1};
    if ( ~isnumeric(n) || n < 0 || round(n) ~= n || ~isscalar(n) )
        error('CHEBFUN:CHEBFUN:dirac:dirac', ...
            'Order of the derivative must be be a non-negative integer.');
    end
    
    if ( n == 0 ) % Trivial case
        d = dirac(f);
        return
    else
        d = diff(dirac(f), n); 
        return
    end
end
    
% Get the epslevel and the domain of f:
tol = epslevel(f);
dom = f.domain;
a = dom(1);
b = dom(end);

% Extract the 'normal' roots of f:
r = roots(f, 'nojump', 'nozerofun');
r = sort(r(:));

% Check roots at the end points of f:
if ( isempty(r) )
    % If there are no roots, still check roots at the end points:
    if ( abs(feval(f, a, 'right')) < 100*tol*f.vscale )
        rootA = 1;
        r = [r; a];
    else
        rootA = 0;
    end
    
    if ( abs(feval(f, b, 'left')) < 100*tol*f.vscale )
        rootB = 1;
        r = [r; b];
    else
        rootB = 0;
    end    
else    
    % If there are roots, check if they are at the end points:
    if ( r(1) > a )
        rootA = 0;
    elseif ( abs(feval(f, a, 'right')) < 100*tol*f.vscale )
        rootA = 1;
        if ( r(1) ~= a )
            r = [a ; r];
        end
    end
    if ( r(end) < b )
        rootB = 0;
    elseif ( abs(feval(f, b, 'left')) < 100*tol*f.vscale )
        rootB = 1;
        if ( r(end) ~= b )
            r = [r ; b];
        end
    end
end

% Initialize a zero CHEBFUN:
d = chebfun(0, [a, b]);

% If there is no root of F within the domain or at the end points, return with a
% zero CHEBFUN:
if ( isempty( r ) )
    return
end

% Check if any of the roots is not simple by looking at the derivative of F:
fp = diff(f);
fpVals = feval(fp, r);
 
% Check root order for interior break-points:
if ( any(abs(fpVals) < 100*tol*fp.vscale) )
    error('CHEBFUN:CHEBFUN:dirac:dirac', ...
        'Function has a root which is not simple');
else
    % Place deltas with appropriate scaling at interior roots.
    deltaMag = 1./abs(fpVals);
end

% Use half of the strength if there is a root at the end point of the
% domain of the input CHEBFUN and update the pointValues:
pointValues = [0; 0];
if ( rootA )
    deltaMag(1) = deltaMag(1)/2;
    pointValues(1) = sign(deltaMag(1))*inf;
end
if ( rootB )
    deltaMag(end) = deltaMag(end)/2;
    pointValues(2) = sign(deltaMag(end))*inf;
end

% Call the DELTAFUN constructor directly:
data.deltaMag = deltaMag.';
data.deltaLoc = r.';
d.funs{1} = deltafun(d.funs{1}, data);
d.pointValues = pointValues;

end
