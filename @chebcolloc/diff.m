function D = diff(disc, m)
%DIFF    Differentiation operator for CHEBCOLLOC discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a COLLOC representation of a Chebyshev polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M (through a better
%   algorithm than multiplication).

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

if ( m == 0 )
    % Trivial case:
    D = eye(sum(disc.dimension));
    
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
    Dn = cell(max(nUnique));
    for nn = nUnique
        Dn{nn} = disc.diffmat(nn, m);
    end
    
    % Scale the blocks appropriately, based on each subinterval length:
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        blocks{k} = Dn{n(k)} * (2/lengths(k))^m;
    end
    
    % Assemble!
    D = blkdiag(blocks{:});
end

end

