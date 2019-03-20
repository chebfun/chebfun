function f = abs( f )
%ABS Absolute value of a BALLFUN.
%   ABS(F) is the absolute value of the BALLFUN F. This function does not work
%   if the function passes through or becomes numerically close to zero.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    return
end 

if ( f.isReal )
    % Positive/negative test.
    % Returns TRUE if there is no sign change.
    
    tol = chebfun2eps;
    bool = false;                  % Assume false

    % Evaluate on a grid:
    X = sample(f);

    X = X(:);

    if ( all( X >= -tol * f.vscale ))   % If all values are nonnegative         
        bool = true;  
    elseif ( all( X <= tol * f.vscale))  % If all values are not positive
        bool = true; 
    end
    if ( ~bool )
        error('CHEBFUN:BALLFUN:abs:notSmooth', ...
            'Sign change detected. Unable to represent the result.'); 
    end
end       

% Still call the constructor in case we missed a change of sign. 
f = compose( f, @abs ); 
end
