function A = vander(f, n)
%VANDER   Vandermonde array-valued CHEBFUN.
%   A = VANDER(F, N) returns a Vandermonde array-valued CHEBFUN whose N columns
%   are powers of the CHEBFUN F, that is A(:,j) = F.^(N-j), j = 0...N-1.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( f(1).isTransposed || numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:vander:row', ...
        'Input must be a scalar-valued column CHEBFUN.')
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
% Reverse the columns to be consistent with MATLAB's built-in vander():
A = fliplr(A);

end
