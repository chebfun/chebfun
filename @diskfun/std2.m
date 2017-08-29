function varargout = std2(varargin)
%STD2   Standard deviation of a DISKFUN.
%   V = STD2(F) computes the standard deviation of a DISKFUN, i.e., 
%
%     STD2(F)^2 = 1/A * sum2(|f(x,y) - m|^2).
%
%   where A is the area of the domain of F.
%
% See also DISKFUN/MEAN, DISKFUN/MEAN2, DISKFUN/STD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = std2@separableApprox(varargin{:});

end