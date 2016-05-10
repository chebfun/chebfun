function G = cheb2quasi(F)
%CHEB2QUASI   Convert an array-valued CHEBFUN to a quasimatrix.
%   CHEB2QUASI(F) converts the array valued CHEBFUN F to a quasimatrix by
%   extracting each of its columns and assigning them to an array of
%   scalar-valued CHEBFUN objects.
%
% See also QUASIMATRIX, QUASI2CHEB, NUM2CELL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(F) > 1 )
    % F is already a quasimatrix!
    G = F;
    return
end

% Collect each column in a cell:
F = num2cell(F);
for k = numel(F):-1:1
    % Assign each column to an element of G:
    G(k) = F{k};
end

end
