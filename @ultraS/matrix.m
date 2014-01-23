function varargout = matrix(disc,dimension,domain)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin > 1 )
    disc.dimension = dimension;
    if ( nargin > 2 )
        disc.domain = domain;
    end
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain)-1) ~= length(disc.dimension) )
    error('Must specify one dimension value for each subinterval.')
end

A = disc.source;
if ( isa(A, 'chebmatrix') )
    c = disc.coeffs;
    outputSpaces = disc.outputSpace;
    L = cell(size(A));
    S = cell(size(A));
    for j = 1:size(A, 1)
        disc.outputSpace = outputSpaces(j);
        for k = 1:size(A, 2)
            disc.coeffs = c{j,k};
            [L{j,k}, S{j,k}] = makeBlocks(disc, A.blocks{j,k});
        end
    end
    if ( isa(A,'linop') )
        [out{1:4}] = applyConstraints(disc,L);
        out{2} = out{2}*cell2mat(S);
    else
        out{1} = cell2mat(L);
    end
    m = max(1,nargout);
    varargout(1:m) = out(1:m);    
else
    disc.coeffs = disc.coeffs{1};
    [varargout{1:nargout}] = makeBlocks(disc, A);
end
end
