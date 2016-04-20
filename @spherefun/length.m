function varargout = length(varargin)
%LENGTH  The rank of a SPHEREFUN.
%   K = LENGTH(F) returns the rank of the SPHEREFUN representation.
%
%   [M, N] = LENGTH( F ) returns the length of the column and 
%   row slices employed in the separable model.  The interpretation 
%   of this quantity depends on the underlying representation. 
%
% See also SPHEREFUN/RANK.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = length@separableApprox(varargin{:});

end
