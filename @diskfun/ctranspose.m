function varargout = ctranspose(varargin)
%'	 Complex conjugate transpose of a DISKFUN.
%   F' is the complex conjugate transpose of F.
%   G = CTRANSPOSE(F) is called for the syntax F'.  
%
%   Since only real-valued DISKFUNS are presently supported, this
%   function is just the same as '.
%
% See also DISKFUN/CONJ, DISKFUN/TRANSPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = ctranspose@separableApprox(varargin{:});

end
