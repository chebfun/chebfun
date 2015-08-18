function varargout = domain(f, varargin)
%DOMAIN   Domain of an ADCHEBFUN.
%   DOMAIN(F) where F is an ADCHEBFUN is the same as DOMAIN(F.FUNC) if F.FUNC is
%   not a scalar, otherwise it returns the domain property of F.
%
% See also CHEBFUN/DOMAIN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isnumeric(f.func) )
    [varargout{1:nargout}] = domain(f.func, varargin{:});
else
    tmp = chebfun(0, f.domain);
    [varargout{1:nargout}] = domain(tmp, varargin{:});
end

end
