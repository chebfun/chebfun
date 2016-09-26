function [A, B] = merge(A, B)
%MERGE   Merge information from two CHEBDSICRETIZATION objects.
%   [A, B] = MERGE(A, B) merges two OPDISCRETIZATIONS A and B.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
% Call the superclass merge:
[A, B] = merge@opDiscretization(A, B);

% Merge the outputSpace:
outputSpace = max(A.outputSpace, B.outputSpace);
A.outputSpace = outputSpace;
B.outputSpace = outputSpace;

end
