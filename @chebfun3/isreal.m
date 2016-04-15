function out = isreal(f)
%ISREAL   Real-valued CHEBFUN3 test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%  
%   ~ISREAL(F) detects CHEBFUN3 object that have an imaginary part even if it is
%   all zero.

if ( isempty(f) )
    out = true;
    return
end

% Get the low rank representation for f. 
[fCore, fCols, fRows, fTubes] = st(f);

% Check individual columns and rows. 
out = isreal(fCore) && isreal(fCols) && isreal(fRows) && isreal(fTubes);

end