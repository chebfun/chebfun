function d = dirac(f, varargin)
% DIRAC delta function
%
% D = DIRAC(F) returns a chebfun D which is zero on the domain of the
% chebfun F except at the simple roots of F, where it is infinite. This
% infinity may be examined by looking at the second row of the matrix
% D.IMPS

% DIRAC(F,N) is the nth derivative of DIRAC(F).

% DIRAC(F) is not defined if F has a zero of order greater than one within
% the domain of F.

% If F has break-points, they should not coinicde with the roots of F.
% However, F can have simple roots at either end points of its domain.

% See also chebfun/heaviside

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty argument:
if ( isempty(f) )
    d = chebfun;
    return
end

% Handle the case for derivatives of delta function
if( nargin > 1 )
    if( nargin > 2 )
        error('CHEBFUN:dirac', 'Too many input arguments');
    end
    % Order of the derivative of dirac delta function
    n = varargin{1};
    if ( ~isa(n, 'double') || n < 0 || round(n) ~= n )
        error('CHEBFUN:dirac', 'order of the derivative must be be a non-negative integer');
    end
    
    if( n == 0 ) % Trivial case
        d = dirac(f);
        return
    else
        d = diff(dirac(f), n); 
        return
    end
end
    
pref = chebpref;
tol = pref.eps;

dom = f.domain;
a = dom(1);
b = dom(end);

% Extract the 'normal' roots of f:
r = roots(f,'nojump','nozerofun');
r = sort(r(:));

% Check roots at the end points of f:
if ( r(1) > a )
    rootA = 0;
elseif ( abs(feval(f, a, 'right')) < 100 * tol* f.vscale )
    rootA = 1;
    if ( r(1) ~= a )
        r = [a; r];
    end
end

if ( r(end) < b )
    rootB = 0;
elseif ( abs(feval(f, b, 'left')) < 100 * tol* f.vscale )
    rootB = 1;
    if( r(end) ~= b )
        r = [r; b];
    end
end

% initialize a zero chebfun
g = chebfun(0, [a, b]);

% If there is no root of F within the domain or at the
% end points, return with a zero chebfun:
if isempty( r )
    return
end

% check if any of the roots is not simple by
% looking at the derivative of F
fp = diff(f);
fpVals = feval(fp, r);
 
% check root order for interior break-points
if ( any(abs(fpVals) < 100 * tol * fp.vscale) )
    error('CHEBFUN:dirac', 'Function has a root which is not simple');
else
    % place delta functions with appropriate scaling at interior roots.
    deltaMag = 1./abs(fpVals);
end

% Use half of the strength if there is a root at the extremal end points of
% the input chebfun:
if ( rootA )
    deltaMag(1) = deltaMag(1)/2;
end

if ( rootB )
    deltaMag(end) = deltaMag(end)/2;
end

% Call the deltafun constructor directly:
d = deltafun(g.funs{1}, deltaMag.', r.');
        
end