function [M, P, B, A, PS] = matrix(disc, dim, domain)
%MATRIX    Convert operator to matrix using COLLOC discretization.
%   MATRIX(DISC) uses the parameters in DISC to discretize DISC.source as a
%   matrix using COLLOC. 
%
%   MATRIX(DISC, DIM, DOMAIN) overrides the native 'dimension' and 'domain'
%   properties in DISC.
%
%   [PA, P, B, A, PS] = MATRIX(...) returns the projection matrix P, the
%   boundary matrix B, the unprojected cell array of square discretizations A,
%   and the projected conversion matrices PS.
%
% See also INSTANTIATE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Parse inputs
if ( nargin > 1 )
    disc.dimension = dim;
end
if ( nargin > 2 )
    disc.domain = domain;
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain) - 1) ~= length(disc.dimension) )
    error('CHEBFUN:chebDiscretisation:matrix:subIntDim', ...
        'Must specify one dimension value for each subinterval.')
end

if ( nargout > 1 && ~isa(disc.source, 'linop') )
    error('CHEBFUN:chebDiscretizartion:matrix', ...
        'MATRIX() of a %s can only return one output.', class(disc.source))
end

% Construct a square representation of each block individually and
% store in a cell array.
[A, S] = instantiate(disc);

% We want output on different format depending on whether the source L is a
% LINOP or something else (typically a standard CHEBMATRIX0:
if ( isa(disc.source, 'linop') )
    
    % Project rows down, and record the projection matrix as well.
    [PA, P, PS] = disc.reduce(A, S);

    % Get constraints:
    B = getConstraints(disc);

    % This should restore squareness to the final matrix.
    M = [ B ; PA ];

else
    
    % Everything should be of the same dimension.
    if ( iscell(A) )
        M = cell2mat(A);
    else 
        M = A;
    end

end


end