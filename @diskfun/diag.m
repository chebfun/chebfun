function F = diag(f,c)
%DIAG(F)   Diagonal of a DISKFUN.
%   G = DIAG(F) returns the CHEBFUN representing g(x) = f(x, x).
%
%   G = diag(F,C) returns the CHEBFUN representing g(x) = f(x, x+c).
%
% See also DISKFUN/TRACE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%[varargout{1:nargout}] = diag@separableApprox(varargin{:});


f.coords = 'polar';

%when c = 0, choose the diagonal radial slice
F = f(pi/4,:); 

%when not(c=0), need to make a chebfun with correct domain: 





end
