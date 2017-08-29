function varargout = std(varargin)
%STD   Standard deviation of a DISKFUN along one variable.
%   G = STD(F) returns the standard deviation of F in the radial variable
%   (default).If F is defined on the rectangle [-pi,pi] x [0,1] then
%
%                         1 
%                        /
%     std(F)^2 = 1/pi    | ( F(THETA,R) - mean(F,1) )^2 dTHETA
%                        /
%                        0
%
%   G = STD(F, FLAG, DIM) takes the standard deviation along the
%   r-variable if DIM = 1 and along the theta-variable (angular) if
%   DIM = 2. The FLAG is ignored and kept in this function so the syntax
%   agrees with the Matlab STD command.
%
% See also CHEBFUN/STD, DISKFUN/MEAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = std@separableApprox(varargin{:});

end