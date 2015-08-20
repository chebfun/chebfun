function C = cat(dim, varargin)
%CAT   Concatenation of CHEBMATRIX objects.
%   C = CAT(dim, A, B,...) concatenates chebmatrices along the indicated
%   dimension following Matlab conventions.
%
% See also CHEBMATRIX.HORZCAT, CHEBMATRIX.VERTCAT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Any singleton operator blocks must be encased in cells to match up with
% CHEBMATRIX blocks.

if ( dim == 1 )
    C = vertcat(varargin{:});
elseif ( dim == 2 )
    C = horzcat(varargin{:});
else
    error('CHEBFUN:CHEBMATRIX:cat:dim', ...
        'Must concatenate along dimension 1 or 2.')
end  

end
