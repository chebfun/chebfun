function varargout = sample(f, varargin)
%SAMPLE   Sample a SINGFUN on an "appropriate" grid.
%   VALUES = SAMPLE(F, N) returns a vector VALUES of the values of F at N
%   "appropriately-chosen" points.  What "appropriately-chosen" means depends
%   on the type of F.SMOOTHPART.
%
%   [VALUES, POINTS] = SAMPLE(F, N) returns also the vector POINTS at which the
%   values were computed.
%
%   [...] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = sample(f.smoothPart, varargin{:});

    % Multiply now with the singular factors:
    if ( f.exponents(1) )
        varargout{1} = varargout{1}.*(1 + x).^f.exponents(1);
    end

    if ( f.exponents(2) )
        varargout{1} = varargout{1}.*(1 - x).^f.exponents(2);
    end
end
