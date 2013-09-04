function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN objects.
%   HORZCAT of a CHEBFUN is not yet supported.

out = chebmatrix(varargin);

% error('CHEBFUN:horzcat:noSupport', 'HORZCAT of a CHEBFUN is not yet supported.');

% [TODO]: Implement this.

end
