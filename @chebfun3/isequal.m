function out = isequal(f, g)
%ISEQUAL   Equality test for CHEBFUN3. 
%   BOOL = ISEQUAL(F, G) returns logical 1 if F and G are same up to relative 
%   machine precision. It returns logical 0 otherwise.
%
% The structure of this code is similar to `dematricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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