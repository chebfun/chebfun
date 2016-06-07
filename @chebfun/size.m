function varargout = size(f, dim)
%SIZE   Size of a CHEBFUN.
%   S = SIZE(F) returns a two-element row vector S = [S1, S2]. If F is a column
%   CHEBFUN, then S1 is infinity and S2 is the number of columns. For a row
%   CHEBFUN, S1 is the number of rows and S2 is infinity. If the finite
%   dimension is > 1, we say F is "array-valued" or a "quasimatrix".
%
%   [S1, S2] = SIZE(F) returns the dimensions of F as separate output variables.
%
%   S = SIZE(F, DIM) returns the dimension specified by the scalar DIM.
%
% See also LENGTH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Infinite dimension:
infDim = inf;

% Number of columns:
numCols = numColumns(f);

% Switch if F is transposed:
if ( f(1).isTransposed )
    tmp = infDim;
    infDim = numCols;
    numCols = tmp;
end

% Parse the output for consistency with MATLAB:
if ( nargout < 2 )
    
    if ( nargin == 1 )
        % Output a vector:
        varargout{1} = [infDim, numCols];
        
    elseif ( dim == 1 )
        % Output first entry:
        varargout{1} = infDim;
        
    elseif ( dim == 2 )
        % Output second entry:
        varargout{1} = numCols;
        
    elseif ( logical(round(dim) - dim) )
        error('CHEBFUN:CHEBFUN:size:dim', ...
            ['Dimension argument must be a positive integer scalar' ...
             ' within indexing range.']);
         
    else
        % We do not allow tensors, so any further dimensions must have size 1.
        varargout{1} = 1;
        
    end
    
else
    
    % Two outputs:
    varargout{1} = infDim;
    varargout{2} = numCols;
    
end

end
