function spy(L, dim, dom, varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

disc = L.discretizer(L);
if ( nargin < 3 )
    dom = L.domain;
end
if ( nargin < 2 )
    dim = 16;
end

numInts = length(dom) - 1;
if ( numel(dim) ~= numInts )
    if ( numel(dim) == 1 )
        dim = repmat(dim, 1, numInts);
    else
        error('Domain/discretisation mismatch.');
    end
end
    
disc.dimension = dim;
disc.domain = dom;
disc = deriveContinuity(disc);
M = matrix(disc);
spy(M, varargin{:});
