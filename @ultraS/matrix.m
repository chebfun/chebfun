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
L = disc.source;
if ( isa(L, 'chebmatrix') )

    % Construct a square representation of each block individually and
    % store in a cell array. The size of the j,k block is determined by
    % disc.dimension + disc.inputDimension(j,k).
    blocks = L.blocks;
    A = cell(size(blocks));    
    for j = 1:size(blocks,1)
        for k = 1:size(blocks,2)
            discJK = extractBlock(disc, j, k);
            A{j,k} = instantiate(discJK);
        end
    end
        
    if ( isa(L, 'linop') )
        [rows, P] = disc.reduce(A);
        PA = cell2mat(rows);
        P = blkdiag(P{:});
        B = getConstraints(disc, A);

        % This should restore squareness to the final matrix.
        M = [ B ; PA ];

    else
        M = cell2mat(A);
    end

        
else
    disc.coeffs = disc.coeffs{1};
    [varargout{1:nargout}] = instantiate(disc, L);
end

end