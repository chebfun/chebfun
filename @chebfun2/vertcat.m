function varargout = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN2 objects.
%
% VERTCAT(F, G) is the vertical concatenation of CHEBFUN2 objects F and G.
% This function returns a CHEBFUN2V object.
%
% [F; G] is equivalent to VERTCAT(F, G).
%
% VERTCAT(F) returns the CHEBFUN2 F.
%
% See also CHEBFUN2V.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = vertcat@separableApprox(varargin{:});

end
