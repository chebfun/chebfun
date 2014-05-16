function G = cell2quasi(F)
%CELL2QUASI   Convert a cell array of CHEBFUN objects to a quasimatrix.
%   CELL2QUASI(F) converts the cell array F of CHEBFUN objects in to a
%   quasimatrix.
%
% See also QUASIMATRIX, CHEB2CELL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~iscell(F) )
    error('CHEBFUN:cell2quasi:notacell', ...
        'Input must be a cell array of CHEBFUN objects.');
end

G = chebfun();
for k = numel(F):-1:1
    G(k) = F{k};
end

end
