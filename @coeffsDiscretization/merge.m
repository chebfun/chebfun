function [varargout] = merge(varargin)
%MERGE   Merge information from two COEFFSDISCRETIZATION objects.
%   [A, B] = MERGE(A, B) merges two COEFFSDISCRETIZATIONS A and B.
%
%   [A1, A2, ...] = MERGE(A1, A2, ...) is the same as above, but for
%   multiple discretizations.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the superclass merge:
[varargout{1:nargout}] = merge@opDiscretization(varargin{:});

% Merge the outputSpace:
outputSpace = 0;
for k = 1:nargin
    outputSpace = max(outputSpace, varargin{k}.outputSpace);
end
for k = 1:nargin
    varargout{k}.outputSpace = outputSpace;
end

end
