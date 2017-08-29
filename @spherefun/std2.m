function varargout = std2(varargin)
%STD2   Standard deviation of a SPHEREFUN.
%   V = STD2(F) computes the standard deviation of a SPHEREFUN, i.e., 
%
%     STD2(F)^2 = 1/(4*pi)*sum2(|F(lambda,theta) - m|^2)
%
%   where m is the mean of F.
%
% See also SPHEREFUN/MEAN, SPHEREFUN/MEAN2, SPHEREFUN/STD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = std2@separableApprox(varargin{:});

end
