function p = pivots(f, str)
%PIVOTS   Pivot values of a CHEBFUN2.
%   PIVOTS(F) returns the pivot values taken during in the constructor by the GE
%   algorithm.
%
%   PIVOTS(F, 'normalize'), returns the normalised pivot values. These numbers
%   are scaled so that the columns and rows have unit 2-norm.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Get the pivot values:
p = f.pivotValues;  

if ( nargin > 1 )
    
    if ( strcmpi(str, 'normalise') || strcmpi(str,'normalize') )
        cscl = norm( f.cols );
        rscl = norm( f.rows );
        % Normalized pivots:
        p = p.*cscl.*rscl;  
    else
        error('CHEBFUN2:PIVOTS:Inputs', 'Unrecognised second argument.');
    end
    
end

% Make column vector:
p = p(:);               

end
