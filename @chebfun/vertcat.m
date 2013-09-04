function out = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN objects.
%   VERTCAT of a CHEBFUN is not yet supported.

out = chebmatrix(varargin.');

% error('CHEBFUN:vertcat:noSupport', 'VERTCAT of a CHEBFUN is not yet supported.');

% [TODO]: Implement this.

end
