function A = vander(f, n)
%VANDER   Vandermonde CHEBFUN quasimatrix.
%   A = VANDER(F, N) returns the Vandermonde quasimatrix whose N columns are
%   powers of the CHEBFUN F, that is A(:,j) = F.^(N-j), j = 0...N-1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( f.isTransposed || min(size(f)) > 1 )
    error('CHEBFUN:vander:row', 'Input must be a scalar-valued column CHEBFUN')
end

% Initialise a cell for storage:
A = cell(1, n);
% The first column (F^0 = 1):
A{1} = chebfun(1, f.domain);
% Loop over the remaining columns:
for j = 1:n-1
    A{j+1} = f.*A{j};
end
% Concatenate the columns to create an array-valued CHEBFUN:
A = horzcat(A{:});

end