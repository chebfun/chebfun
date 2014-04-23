function [M, P, B, A] = matrix(disc, dim, domain)
%MATRIX    Convert operator to matrix using COLLOC discretization.
%   MATRIX(DISC) uses the parameters in DISC to discretize DISC.source as a
%   matrix using COLLOC. 
%
%   MATRIX(DISC, DIM, DOMAIN) overrides the native 'dimension' and 'domain'
%   properties in DISC.
%
%   [PA, P, B, A] = MATRIX(...) returns the additional component matrices
%   resulting from boundary condition manipulations, as described in
%   APPLYCONSTRAINTS.
%
%   See also: INSTANTIATE

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Parse inputs
if ( nargin > 1 )
    disc.dimension = dim;
    if ( nargin > 2 )
        disc.domain = domain;
    end
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain) - 1) ~= length(disc.dimension) )
    error('Must specify one dimension value for each subinterval.')
end

L = disc.source;
if ( isa(L, 'chebmatrix') )
    
    % TODO: Do these things ever occur?
    if ( isempty(disc.inputDimension) )
        disc.inputDimension = 0;
    end
    if ( numel(disc.inputDimension) == 1 )
        disc.inputDimension = repmat(disc.inputDimension, size(L.blocks));
    end
    
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

    % We want output on different format depending on whether the source L is a
    % LINOP or another object (most likely a CHEBMATRIX).
    if ( isa(L, 'linop') )
        % Project rows down, and record the projection matrix as well.
        [rows, P] = disc.reduce(A);
        PA = cell2mat(rows);
        P = blkdiag(P{:});
        B = getConstraints(disc, L.blocks);
        
        % This should restore squareness to the final matrix.
        M = [ B; PA ];
    else
        M = cell2mat(A);
    end

else
    % The source must be a chebfun, or...?
    % Note, this is called by ultraS for functionalBlocks
    [varargout{1:nargout}] = instantiate(disc,L);
    M = instantiate(disc,L);
end

end