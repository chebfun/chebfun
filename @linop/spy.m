function spy(L, dim, dom, varargin)

disc = L.discretizer(L);
if ( nargin < 3 )
    dom = L.domain;
end
if ( nargin < 2 )
    dim = 100;
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
