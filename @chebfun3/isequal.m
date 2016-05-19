function out = isequal(f, g)
%ISEQUAL   Equality test for CHEBFUN3. 
%   BOOL = ISEQUAL(F, G) returns 0 or 1. If returns 1 then F and G are the 
%   same CHEBFUN3 objects, up to relative machine precision. It returns 0 
%   if F and G are not the same up to relative machine precision. 

% The structure of this code is similar to `dematricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    if ( isempty(g) )
        out = true; 
    else
        out = false; 
    end
    return
end

% Get low rank representation of f: 
[fCore, fCols, fRows, fTubes] = tucker(f);

% Get low rank representation of g:
[gCore, gCols, gRows, gTubes] = tucker(g);

% Test every part: 
out = ( isequal(fCore, gCore) & isequal(fCols, gCols) & ...
    isequal(fRows, gRows) & isequal(fTubes, gTubes) );

end