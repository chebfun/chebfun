function varargout = mean2(varargin)
%MEAN2   Mean of a CHEBFUN2.
%   V = MEAN2(F) returns the mean of a CHEBFUN2:
%
%   V = 1/A*integral2( f )
%
% 	where the A is the area of the domain of F.
%
% See also MEAN, STD2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mean2@separableApprox(varargin{:});

end
