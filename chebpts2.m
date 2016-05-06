function varargout = chebpts2(varargin)
%CHEBPTS2   Chebyshev tensor product grid.
%   [XX, YY] = CHEBPTS2(N) constructs an N by N grid of Chebyshev tensor points
%   on [-1 1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY) constructs an NX by NY grid of Chebyshev tensor
%   points on [-1 1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY, [a b c d]) constructs an NX by NY grid of
%   Chebyshev tensor points on the rectangle [a b] x [c d].
%
%   [XX, YY] = CHEBPTS2(NX, NY, D, KIND) constructor Chebyshev tensor grid of
%   the kind KIND. KIND = 2 is default. 
% 
% See also CHEBPTS, CHEBFUN2.CHEBPTS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This method is just a wrapper for CHEBFUN2.CHEBPTS2().
[varargout{1:nargout}] = chebfun2.chebpts2(varargin{:}); 

end 
