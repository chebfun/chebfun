function X = mrdivide(B, A)
%/	Right matrix divide for a FUNCHEB2.
%
%   B/A divides the FUNCHEB2 B by a scalar A. More generally, it gives the least
%   squares solution (wrt the continuous L^2 norm) to X*A = B when either A or B
%   is a FUNCHEB2. NOTE: In the latter case, formally it is X.' that is returned
%   (as FUNCHEB2 objects are always columns).
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case A is a matrix: X*A = B  ==> X = B/A = (Q*R)/A = Q*(R/A) 
%
% Case B is a matrix: X*A = X*(Q*R) = B ==> X = (B/R)*Q' ==> X' = Q*(B/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( size(A, 2) ~= size(B, 2) && ~( isa(A, 'double') && isscalar(A) ) )
    
    error('CHEBFUN:FUNCHEB2:mrdivide:size', 'Matrix dimensions must agree.');
            
elseif ( isa(A, 'double') )  % funcheb2 / double

    if ( ~any(A(:)) )  
        X = funcheb2(NaN(1, size(B, 2)));
    elseif ( isscalar(A) )
        % Scalar case is easy:
        X = B;                              % Copy B to X
        X.values = B.values/A;              % Divide values
        X.coeffs = B.coeffs/A;              % Divide coeffs
        X.vscale = B.vscale/abs(A);         % Divide vscale
    else
        % For matrix case, we do least squares via QR:
        [Q, R] = qr(B, 0);
        X = Q*(R/A);
    end
    
        
elseif ( isa(B, 'double') )  % double / funcheb2
    
    % Here B is a double and A is a FUNCHEB2. Do least squares via QR:
    [Q, R] = qr(A, 0);
    
    % Return the transpose for the output.
    X = Q*(B/R).';
    
elseif ( isa(A, 'funcheb2') && isa(B, 'funcheb2') )
    
    error('CHEBFUN:FUNCHEB2:mrdivide:funfun', ...
        'Use ./ to divide by a FUNCHEB2.');

else
    
    error('CHEBFUN:FUNCHEB2:mrdivide:derp', '%s/%s is not well-defined.', ...
        class(B), class(A));

end