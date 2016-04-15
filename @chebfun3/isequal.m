function out = isequal(f, g)
%ISEQUAL    Equality test for CHEBFUN3.
% 
%   BOOL = ISEQUAL(F, G) returns 0 or 1. If returns 1 then F and G are the 
%   same CHEBFUN3 objects, up to relative machine precision. If returns 0 
%   then F and G are not the same up to relative machine precision. 

% The structure of this code is similar to `dematricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

if ( isempty(f) )
    if ( isempty(g) )
        out = true; 
    else
        out = false; 
    end
    return
end

% Get the low rank representation for f. 
[fCore, fCols, fRows, fTubes] = st(f);

% Get the low rank representation for g. 
[gCore, gCols, gRows, gTubes] = st(g);

% Test every part: 
out = ( isequal(fCore, gCore) & isequal(fCols, gCols) & ...
    isequal(fRows, gRows) & isequal(fTubes, gTubes) );

end