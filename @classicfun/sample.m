function varargout = sample(f, varargin)
%SAMPLE   Sample a CLASSICFUN on an "appropriate" grid.
%   VALUES = SAMPLE(F, N) returns a vector VALUES of the values of F at N
%   "appropriately-chosen" points.  What "appropriately-chosen" means depends
%   on the type of representation on which F is ultimately based.
%
%   [VALUES, POINTS] = SAMPLE(F, N) returns also the vector POINTS at which the
%   values were computed.
%
%   [...] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = sample(f.onefun, varargin{:});
    if ( nargout > 1 )
        varargout{2} = f.mapping.For(varargout{2});
    end
end
