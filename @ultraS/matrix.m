function L = matrix(disc,dimension,domain)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% TODO: error checking on inputs
if ( nargin > 1 )
    disc.dimension = dimension;
    if ( nargin > 2 )
        disc.domain = domain;
    end
end

A = disc.source;
%            validate(disc);
if ( isa(A, 'chebmatrix') )
    c = disc.coeffs;
    outputSpaces = disc.outputSpace;
    L = cell(size(A));
    for j = 1:size(A, 1)
        disc.outputSpace = outputSpaces(j);
        for k = 1:size(A, 2)
            disc.coeffs = c{j,k};
            L{j,k} = blockDiscretize(disc, A.blocks{j,k});
        end
    end
    if ( isa(A,'linop') )
        L = useConstraints(disc,L);
    end
else
    disc.coeffs = disc.coeffs{1};
    L = blockDiscretize(disc, A);
end
end
