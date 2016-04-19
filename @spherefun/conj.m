function varargout = conj(varargin)
%CONJ   Complex conjugate of a SPHEREFUN.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
%   Since only real-valued SPHEREFUNS are presently supported, this
%   function is trivial.
%
% See also SPHEREFUN/REAL, SPHEREFUN/IMAG. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = conj@separableApprox(varargin{:});

end
