function f = extractColumns(f, colIdx)
%EXTRACTCOLUMNS   Extract columns (or rows) of SINGFUN.
%   EXTRACTCOLUMNS(F, COLIDX) returns the SINGFUN F when COLIDX = 1 and throws
%   an error otherwise (as SINGFUNs may not be array-valued).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(colIdx) ~= 1 || colIdx ~= 1 )
    error('CHEBFUN:SINGFUN:extractColumns:dim', ...
        'Index exceeds matrix dimensions.');
end   

% DEVELOPER NOTE: Singfuns currently only have a single column.
    
end
