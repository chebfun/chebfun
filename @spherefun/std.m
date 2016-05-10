function varargout = std(varargin)
%STD   Standard deviation of a SPHEREFUN along one variable.
%   G = STD(F) returns the standard deviation of F in the theta-variable
%   (or latitude), which is the default.
%   That is, if F is defined on the rectangle [-pi,pi] x [0,pi] then
%
%                         pi 
%                        /
%     std(F)^2 = 1/pi    | ( F(LAMBDA,THETA) - mean(F,1) )^2 dTHETA
%                        /
%                        0
%
%   G = STD(F, FLAG, DIM) takes the standard deviation along the
%   theta-variable if DIM = 1 and along the lambda-variable (longitude) if
%   DIM = 2. The FLAG is ignored and kept in this function so the syntax
%   agrees with the Matlab STD command.
%
% See also CHEBFUN/STD, SPHEREFUN/MEAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = std@separableApprox(varargin{:});

end
