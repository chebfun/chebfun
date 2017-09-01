function varargout = mean2(varargin)
%MEAN2   Mean of a SPHEREFUN.
%   V = MEAN2(F) returns the mean of a SPHEREFUN: 
% 
%   V = 1/A*integral2( f )
% 
% 	where the A is the area of the sphere, i.e. 4*pi. 
%
% See also SPHEREFUN/MEAN, SPHEREFUN/STD2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mean2@separableApprox(varargin{:});

end
