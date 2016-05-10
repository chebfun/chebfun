function varargout = sum(varargin)
%SUM   Definite Integration of a CHEBFUN2.
%   G = sum(F,DIM) where DIM is 1 or 2 integrates only over Y or X respectively,
%   and returns as its output a chebfun in the remaining variable.
%
%   G = sum(F) is the same as sum(F,1)
%
% See also SUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sum@separableApprox(varargin{:});

end
