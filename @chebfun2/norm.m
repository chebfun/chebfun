function varargout = norm(varargin)
%NORM       Norm of a CHEBFUN2 object.
% For CHEBFUN2 objects:
%    NORM(F) = sqrt(integral of abs(F)^2).
%    NORM(F, 2) = largest singular value of F.
%    NORM(F,'fro') is the same as NORM(F).
%    NORM(F, 1) = NOT IMPLEMENTED.
%    NORM(F, inf) = global maximum in absolute value.
%    NORM(F, max) = global maximum in absolute value.
%    NORM(F, min) = NOT IMPLEMENTED
%
% Furthermore, the inf norm for CHEBFUN2 objects also returns a second output,
% giving a position where the max occurs.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = norm@separableApprox(varargin{:});

end
