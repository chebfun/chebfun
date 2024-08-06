function D = diff(disc, m)
%DIFF    Differentiation operator for CHEBCOLLOC discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a COLLOC representation of a Chebyshev polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M (through a better
%   algorithm than multiplication).

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

pref = cheboppref;

if ( m == 0 )
    % Trivial case:
    if ( pref.sparse )
        D = speye(sum(disc.dimension));
    else
        D = eye(sum(disc.dimension));
    end

else
    % Find the diagonal blocks.
    
    % Store information about domain and dimensions.
    numIntervals = disc.numIntervals;   
    n = disc.dimension;
    d = disc.domain;
    lengths = diff(d);
    
    % Check how many unique discretization sizes we've got:
    nUnique = unique(n);
    
    % Create a discretization matrix for each unique discretization size:
    Dn = cell(max(nUnique),1);
    for nn = nUnique
        Dn{nn} = disc.diffmat(nn, m);
    end
    if ( pref.sparse ) % Sparsify
        for nn = nUnique
            Dn{nn} = sparse(Dn{nn});
        end
    end

    % Scale the blocks appropriately, based on each subinterval length:
    blocks = cell(numIntervals,1);
    for k = 1:numIntervals
        blocks{k} = Dn{n(k)} * (2/lengths(k))^m;
    end

    % Assemble!
    % D = blkdiag(blocks{:}); % <-- Lots of overhead here.
    D = matlab.internal.math.blkdiag(blocks{:});
end

end

