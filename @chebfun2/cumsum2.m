function varargout = cumsum2(varargin)
%CUMSUM2   Double indefinite integral of a CHEBFUN2.
%   F = CUMSUM2(F) returns the double indefinite integral of a CHEBFUN2. That is
%                   y  x
%                  /  /
%   CUMSUM2(F) =  |  |   f(x,y) dx dy   for  (x,y) in [a,b] x [c,d],
%                 /  /
%                c  a
%
%   where [a,b] x [c,d] is the domain of f.
%
% See also CUMSUM, SUM, SUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cumsum2@separableApprox(varargin{:});

end
