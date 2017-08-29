function varargout = ctranspose(varargin)
%'	 Complex conjugate transpose of a SPHEREFUN.
%   F' is the complex conjugate transpose of F.
%   G = CTRANSPOSE(F) is called for the syntax F'.  
%
%   Since only real-valued SPHEREFUNS are presently supported, this
%   function is just the same as '.
%
% See also SPHEREFUN/CONJ, SPHEREFUN/TRANSPOSE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = ctranspose@separableApprox(varargin{:});

end
