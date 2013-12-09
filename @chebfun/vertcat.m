function out = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN objects.
%   VERTCAT of a CHEBFUN is not yet supported.

% TODO: Document.
% TODO: Test.

out = chebmatrix(varargin.');

end
