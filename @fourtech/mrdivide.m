function X = mrdivide(A, B)
%/   Right matrix divide for a FOURTECH.
%
%   A/B divides the FOURTECH A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is a FOURTECH.  Note that in the latter case, formally
%   it is X.' that is returned, as FOURTECH objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case B is a matrix: X*B = A  ==> X = A/B = (Q*R)/B = Q*(R/A) 
%
% Case A is a matrix: X*B = X*(Q*R) = A ==> X = (A/R)*Q' ==> X' = Q*(A/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(B, 2) ~= size(A, 2)) && ~(isa(B, 'double') && isscalar(B)) )
    error('CHEBFUN:FOURTECH:mrdivide:size', 'Matrix dimensions must agree.');
elseif ( isa(B, 'double') )  % FOURTECH / double
    if ( ~any(B(:)) )  
        X = A.make(NaN(1, size(A, 2)));
    elseif ( isscalar(B) )
        % Scalar case is easy:
        X = A;                              % Copy A to X
        X.values = A.values/B;              % Divide values
        X.coeffs = A.coeffs/B;              % Divide coeffs
        X.vscale = A.vscale/abs(B);         % Divide vscale
    else
%         % For matrix case, we do least squares via QR:
%         [Q, R] = qr(A, 0);
%         X = Q*(R/B);
        error('CHEBFUN:FOURTECH:mrdivide','This option requires the QR which is not implemented for FOURTECH objects yet');
    end
elseif ( isa(A, 'double') )  % double / FOURTECH
    error('CHEBFUN:FOURTECH:mrdivide','This option requires the QR which is not implemented for FOURTECH objects yet');
%     % Here A is a double and B is a FOURTECH. Do least squares via QR:
%     [Q, R] = qr(B, 0);
%     
%     % Return the transpose for the output.
%     X = Q*(A/R).';
elseif ( isa(B, 'fourtech') && isa(A, 'fourtech') )
    error('CHEBFUN:FOURTECH:mrdivide:chebtechDivChebtech', ...
        'Use ./ to divide by a FOURTECH.');
else
    error('CHEBFUN:FOURTECH:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));
end

end
