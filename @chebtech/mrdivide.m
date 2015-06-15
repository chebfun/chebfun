function X = mrdivide(A, B)
%/   Right matrix divide for a CHEBTECH.
%
%   A/B divides the CHEBTECH A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is a CHEBTECH.  Note that in the latter case, formally
%   it is X.' that is returned, as CHEBTECH objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case B is a matrix: X*B = A  ==> X = A/B = (Q*R)/B = Q*(R/A) 
%
% Case A is a matrix: X*B = X*(Q*R) = A ==> X = (A/R)*Q' ==> X' = Q*(A/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(B, 2) ~= size(A, 2)) && ~(isa(B, 'double') && isscalar(B)) )
    error('CHEBFUN:CHEBTECH:mrdivide:size', 'Matrix dimensions must agree.');
elseif ( isa(B, 'double') )  % CHEBTECH / double
    if ( ~any(B(:)) )  
        X = A.make(NaN(1, size(A, 2)));
    elseif ( isscalar(B) )
        % Scalar case is easy:
        X = A;                              % Copy A to X
        X.coeffs = A.coeffs/B;              % Divide coeffs
    else
        % For matrix case, we do least squares via QR:
        [Q, R] = qr(A, 0);
        X = Q*(R/B);
    end
elseif ( isa(A, 'double') )  % double / CHEBTECH
    % Here A is a double and B is a CHEBTECH. Do least squares via QR:
    [Q, R] = qr(B, 0);
    
    % Return the transpose for the output.
    X = Q*(A/R).';
elseif ( isa(B, 'chebtech') && isa(A, 'chebtech') )
    error('CHEBFUN:CHEBTECH:mrdivide:chebtechDivChebtech', ...
        'Use ./ to divide by a CHEBTECH.');
else
    error('CHEBFUN:CHEBTECH:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));

end

end
