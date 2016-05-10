function X = mrdivide(A, B)
%/   Right matrix divide for a TRIGTECH.
%   A/B divides the TRIGTECH A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is a TRIGTECH.  Note that in the latter case, formally
%   it is X.' that is returned, as TRIGTECH objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case B is a matrix: X*B = A  ==> X = A/B = (Q*R)/B = Q*(R/A) 
%
% Case A is a matrix: X*B = X*(Q*R) = A ==> X = (A/R)*Q' ==> X' = Q*(A/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(B, 2) ~= size(A, 2)) && ~(isa(B, 'double') && isscalar(B)) )
    error('CHEBFUN:TRIGTECH:mrdivide:size', 'Matrix dimensions must agree.');
    
elseif ( isa(B, 'double') )  % TRIGTECH / double.
    if ( ~any(B(:)) )  
        X = A.make(NaN(1, size(A, 2)));
    elseif ( isscalar(B) )
        % Scalar case is easy:
        X = A;                              % Copy A to X.
        X.values = A.values/B;              % Divide values.
        X.coeffs = A.coeffs/B;              % Divide coeffs.
        X.isReal = A.isReal & isreal(B);
    else
        % For matrix case, we do least squares via QR:
        [Q, R] = qr(A, 0);
        X = Q*(R/B);
    end
    
elseif ( isa(A, 'double') )  % double / TRIGTECH.
    % Here A is a double and B is a TRIGTECH. Do least squares via QR:
    [Q, R] = qr(B, 0);
    
    % Return the transpose for the output.
    X = Q*(A/R).';
    
elseif ( isa(B, 'trigtech') && isa(A, 'trigtech') )
    error('CHEBFUN:TRIGTECH:mrdivide:trigtechDivTrigtech', ...
        'Use ./ to divide by a TRIGTECH.');
else
    error('CHEBFUN:TRIGTECH:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));
end

end
