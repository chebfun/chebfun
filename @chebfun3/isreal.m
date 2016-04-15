function out = isreal(f)
%ISREAL   Real-valued CHEBFUN3 test.
%   ISREAL(F) returns logical true if F does not have an nonzero imaginary
%   part and false otherwise.
%  
%   (This is slightly different from the Matlab convention, where
%   isreal(x) is false if x is a complex number whose imaginary part is 0.)

if ( isempty(f) )
    out = true;
    return
end

% Get the low rank representation for f. 
[fCore, fCols, fRows, fTubes] = st(f);

% Check individual columns and rows. 
out = isreal(fCore) && isreal(fCols) && isreal(fRows) && isreal(fTubes);

end
