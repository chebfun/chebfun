function varargout = matrix(disc, dimension, domain)
%MATRIX    Convert operator to matrix using ULTRAS discretization.
%   MATRIX(DISC) uses the parameters in DISC to discretize DISC.source as a
%   matrix. 
%
%   MATRIX(DISC, DIM, DOMAIN) overrides the native 'dimension' and 'domain'
%   properties in DISC.
%
%   [PA, P, B, A] = MATRIX(...) returns the additional component matrices
%   resulting from boundary condition manipulations, as described in
%   ULTRAS.APPLYCONSTRAINTS.
%
%   See also: ULTRAS/INSTANTIATE.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin > 1 )
    disc.dimension = dimension;
    if ( nargin > 2 )
        disc.domain = domain;
    end
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain) - 1) ~= length(disc.dimension) )
    error('Must specify one dimension value for each subinterval.')
end

% TODO: Add some documentation below please. AB, 1/3/14.
A = disc.source;
if ( isa(A, 'chebmatrix') )
    c = disc.coeffs;
    outputSpaces = disc.outputSpace;
    L = cell(size(A));
    S = cell(size(A));
    for j = 1:size(A, 1)
        disc.outputSpace = outputSpaces(j);
        for k = 1:size(A, 2)
            disc.coeffs = c{j, k};
            [L{j, k}, S{j, k}] = instantiate(disc, A.blocks{j, k});
        end
    end
    if ( isa(A, 'linop') )
        [out{1:4}] = applyConstraints(disc, L);
        out{2} = out{2}*cell2mat(S);
    else
        out{1} = cell2mat(L);
    end
    m = max(1, nargout);
    varargout(1:m) = out(1:m);    
else
    disc.coeffs = disc.coeffs{1};
    [varargout{1:nargout}] = instantiate(disc, A);
end
end
