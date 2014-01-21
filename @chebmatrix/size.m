function varargout = size(L, varargin)
%SIZE      Number of blocks in each dimension.
%   S = SIZE(L) returns both dimensions.
%   S = SIZE(L, K) returns Kth dimension (K = 1,2).
%   [M, N] = SIZE(L) returns both as scalars.
% 
%   See also CHEBMATRIX.BLOCKSIZES.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

[varargout{1:nargout}] = size(L.blocks, varargin{:});

end
