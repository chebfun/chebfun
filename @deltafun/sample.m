function [values, points] = sample(f, varargin)
%SAMPLE   Sample a DELTAFUN on an "appropriate" grid.
%   VALUES = SAMPLE(F, N) returns a vector VALUES of the values of F at N
%   "appropriately-chosen" points.  What "appropriately-chosen" means depends
%   on the type of SMOOTHFUN object on which F is ultimately based.
%
%   [VALUES, POINTS] = SAMPLE(F, N) returns also the vector POINTS at which the
%   values were computed.
%
%   [...] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = sample(f.funPart, varargin{:});
end
