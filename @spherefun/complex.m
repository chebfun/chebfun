function varargout = complex(varargin)
%COMPLEX Construct complex SPHEREFUN2 from real and imaginary parts.
%   C = COMPLEX(A, B) returns the complex SPHEREFUN2 A + Bi, where A and B are
%   real valued SPHEREFUN2 objects with the same domain.
%
%   C = COMPLEX(A) for real SPHEREFUN2 A returns the complex result C with real
%   part A and all zero imaginary part. isreal(C) returns false.
%
% See also IMAG, CONJ, ABS, REAL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = complex@separableApprox(varargin{:});
end
