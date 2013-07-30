function f = cumsum(f, m, pref)
% CUMSUM   Indefinite integral.
%
% G = CUMSUM(F) is the indefinite integral of the chebfun F. G will typically be
% normalised so that G(F.domain(1)) = 0. The exception to this is when computing
% indefinite integrals of functions with exponents less than minus 1. In this
% case, the arbitrary constant in the indefinite integral is chosen to make the
% representation of G as simple as possible. Dirac deltas already existing in F
% will decrease their degree.
%
% CUMSUM(F, M) returns the Mth integral of F. If N is not an integer CUMSUM(F,
% N) returns the fractional integral of order N as defined by the
% Riemann-Liouville integral.
%
% CUMSUM does not currently support chebfuns whose indefinite integral diverges
% (i.e. has exponents <-1) when using nontrivial maps. Even for chebfuns with a
% bounded definite integral, nontrivial maps will be slow.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% [TODO]: m > 1.

% Trivial case:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    m = 1;
    pref = chebfun.pref;
elseif ( nargin == 2 )
    if ( isstruct(m) )
        pref = m;
        m = 1;
    else
        pref = chebfun.pref;
    end
end

dom = f.domain;
funs = f.funs;
nfuns = numel(funs);

if ( size(f.impulses, 1) > 1 )
    imps = f.impulses(2,:);
else
    imps = zeros(size(dom));
end

fa = imps(1);
for j = 1:nfuns
    csfj = cumsum(funs{j});

%     if ( nfuns > 1 )
%     % This is because unbounded functions may not be zero at left.
%         lval = get(csfj, 'lval');
%         if ( ~isinf(lval) && ~isnan(lval) )
%             csfj = csfj - lval;
%         end
%     end

    funs{j} = csfj + fa;
    fa = get(funs{j},'rval') + imps(j+1);
end

newimps = zeros(1, nfuns+1);
for j = 1:nfuns
    newimps(j) = get(funs{j}, 'lval');
end
newimps(nfuns+1) = get(funs{nfuns}, 'rval');

f.funs = funs;
f.impulses = [ newimps ; f.impulses(3:end,:) ];

end


