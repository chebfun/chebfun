function varargout = conj(varargin)
%CONJ   Complex conjugate of a DISKFUN.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
%   Since only real-valued DISKFUNS are presently supported, this
%   function is trivial.
%
% See also DISKFUN/REAL, DISKFUN/IMAG. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = conj@separableApprox(varargin{:});

end
