function varargout = std(varargin)
%STD   Standard deviation of a CHEBFUN2 along one variable.
%   G = STD(F) returns the standard deviation of F in the y-variable (default).
%   That is, if F is defined on the rectangle [a,b] x [c,d] then
%
%                         d
%                        /
%     std(F)^2 = 1/(d-c) | ( F(x,y) - mean(F,1) )^2 dy
%                        /
%                        c
%
%   G = STD(F, FLAG, DIM) takes the standard deviation along the y-variable if
%   DIM = 1 and along the x-variable if DIM = 2. FLAG is ignored and kept in
%   this function so the syntax agrees with that of the Matlab STD command.
%
% See also CHEBFUN/STD, CHEBFUN2/MEAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = std@separableApprox(varargin{:});

end
