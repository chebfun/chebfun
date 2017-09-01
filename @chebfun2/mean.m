function varargout = mean(varargin)
%MEAN   Average or mean value of a CHEBFUN2.
%   MEAN(F) takes the mean in the y-direction (default), i.e.,
%          MEAN(F) = 1/(ymax-ymin) sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%   y-direction and if DIM = 2 then it is the x-direction.
%
% See also MEAN2, STD2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mean@separableApprox(varargin{:});

end
