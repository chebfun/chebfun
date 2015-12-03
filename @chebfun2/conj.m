function varargout = conj(varargin)
%CONJ   Complex conjugate of a CHEBFUN2.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = conj@separableApprox(varargin{:});

end
