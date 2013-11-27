function G = repmat(f, M, N)
%REPMAT   Replicate and tile a CHEBFUN.
%   REPMAT(F, M, N) or REPMAT(F, [M, N]) creates a CHEBFUN quasimatrix by tiling
%   copies of F. If F is a column quasimatrix, then REPMAT(F, 1, N) returns a
%   quasimatrix with N*size(F, 2) chebfun columns. If F is a row quasimatrix,
%   REPMAT(F, M, 1) returns a quasimatrix with M*size(F, 1).

% See also HORZCAT, VERTCAT, CAT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin == 2 )
    if ( length(M) ~= 2 )
        error('CHEBFUN:repmat:NotEnoughInputs', ...
            'Requires REPMAT(F, M, N) or REPMAT(F, [M, N]).')
    end
    N = M(2);  
    M = M(1);
end

if ( f.isTransposed )
    % REPMAT a row CHEBFUN:
    if ( N ~= 1 )
        error('CHEBFUN:repmat:row',...
            'Use REPMAT(F, M, 1) to replicate and tile row CHEBFUN objects.')
    else
        G = cell(1,M);
        for j = 1:M, 
            G{j} = f;
        end
        G = vertcat(G{:});
    end
else
    % REPMAT a column CHEBFUN:
    if ( M ~= 1 )
        error('CHEBFUNchebfun:repmat:col',...
            'Use REPMAT(F, 1, N) to replicate and tile column CHEBFUN objects.')
    else
        G = cell(1,N);
        for j = 1:N
            G{j} = f;
        end
        G = horzcat(G{:});
    end
end