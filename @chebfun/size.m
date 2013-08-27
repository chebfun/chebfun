function varargout = size(f, dim)
%SIZE   Size of a CHEBFUN.
%   [S1, S2] = SIZE(F) returns S1, the number of piecewise smooth components of
%   F, and S2, the number of columns in F. If S2 > 1, we say that F is
%   "array-valued".
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Number of FUNs:
numFuns = numel(f.funs);

% Number of columns:
if ( numFuns ~= 0 )
    numCols = size(f.funs{1}, 2);
else
    numCols = 0;
end

% Switch if F is transposed:
if ( f.isTransposed )
    tmp = numFuns;
    numFuns = numCols;
    numCols = tmp;
end

% Parse the output for consistency with MATLAB:
if ( nargout < 2 )
    if ( nargin == 1 )
        % Output a vector:
        varargout{1} = [numFuns, numCols];
    elseif ( dim == 1 )
        % Output first entry:
        varargout{1} = numFuns;
    elseif ( dim == 2 )
        % Output second entry:
        varargout{1} = numCols;
    elseif ( logical(round(dim) - dim) )
        error('CHEBFUN:size:dim', ...
            ['Dimension argument must be a positive integer scalar' ...
             ' within indexing range.']);
    else
        % We do not allow tensors, so any further dimensions must have size 1.
        varargout{1} = 1;
    end
else
    % Two outputs:
    varargout{1} = numFuns;
    varargout{2} = numCols;
end

end
