function out = isequal( f, g )
%ISEQUAL Equality test for SEPARABLEAPPROX.  
% 
% BOL = ISEQUAL(F, G) returns 0 or 1. If returns 1 then F and G are the same
% SEPARABLEAPPROX, up to relative machine precision. If returns 0 then F and G are
% not the same up to relative machine precision. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    if ( isempty( g ) )
        out = true; 
    else
        out = false; 
    end
    return
end

% Get the low rank representation for f. 
fcols = f.cols; 
frows = f.rows; 
fpiv = f.pivotValues;

% Get the low rank representation for g. 
gcols = g.cols; 
grows = g.rows; 
gpiv = g.pivotValues;

% Test every part: 
out = ( isequal(fcols, gcols) & isequal(frows, grows) & isequal(fpiv, gpiv) );

end
