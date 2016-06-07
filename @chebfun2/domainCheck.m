function varargout = domainCheck(varargin)
%DOMAINCHECK   True if the domains of two CHEBFUN2 objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the domains of the two
%   CHEBFUN2 objects F and G coincide up to a tolerance depending on their
%   horizontal scales or if both F and G are empty CHEBFUN objects.
%
% See also CHEBFUN/DOMAINCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = domainCheck@separableApprox(varargin{:});

end
