function bol = isequal( f, g )
%ISEQUAL Equality test for chebfun2.  
% 
% BOL = ISEQUAL(F,G) returns 0 or 1. If returns 1 then F and G are the same
% chebfun2, up to relative machine precision. If returns 0 then F and G are
% not the same up to relative machine precision. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) )
    if ( isempty(g) )
        bol = 1; 
    else
        bol = 0; 
    end
    return
end

% Get the low rank representation for f. 
fcols = f.cols; 
frows = f.rows; 
fpiv = f.pivotValues;

% Get the low rank representation for g. 
gcols = f.cols; 
grows = f.rows; 
gpiv = f.pivotValues;

% Test: 
bol = ( isequal(fcols, gcols) & isequal(frows, grows) & isequal(fpiv, gpiv) );

end