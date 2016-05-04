function out = numColumns(f)
%NUMCOLUMNS   Number of columns (or rows) of a CHEBFUN quasimatrix.
%   NUMCOLUMNS(F) returns the number of columns of a column CHEBFUN, or the
%   number of rows of a row CHEBFUN. It is equivalent to min(size(f)).
%
% See also SIZE, NUMEL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    % An empty CHENFUN has no columns:
    out = 0;
elseif ( numel(f) == 1 )
    % Possible array-valued:
    out = size(f.funs{1}, 2);
else
    % A quasimatrix:
    out = numel(f);
end

end
