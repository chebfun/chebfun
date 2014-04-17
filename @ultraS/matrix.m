function [M, P, B, A] = matrix(disc, dimension, domain)
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
    inputDimension = disc.inputDimension;
    L = cell(size(A));
    S = cell(size(A));
    for j = 1:size(A, 1)
        disc.outputSpace = outputSpaces(j);
        for k = 1:size(A, 2)
            disc.coeffs = c{j, k};
            disc.inputDimension = inputDimension(j,k);
            [L{j,k}, S{j,k}] = instantiate(disc, A.blocks{j, k});
        end
    end
    disc.inputDimension = inputDimension;
    
    if ( isa(A, 'linop') )
        [rows, P] = disc.reduce(L);
        PA = cell2mat(rows);
        P = blkdiag(P{:});
        B = getConstraints(disc, L);

        % This should restore squareness to the final matrix.
        M = [ B ; PA ];

        if ( size(M, 1) ~= size(M, 2) )
            % TODO: Improve this warning.
             warning('Matrix is not square!');
        end
    else
        M = cell2mat(L);
    end

    

    

        
else
    disc.coeffs = disc.coeffs{1};
    [varargout{1:nargout}] = instantiate(disc, A);
end

end
