function X = mrdivide(B, A)
%/   Right matrix divide for a CHEBTECH.
%
%   B/A divides the CHEBTECH B by a scalar A. More generally, it gives the least-
%   squares solution (with respect to the continuous L^2 norm) to X*A = B when
%   either A or B is a CHEBTECH.  Note that in the latter case, formally it is
%   X.' that is returned, as CHEBTECH objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case A is a matrix: X*A = B  ==> X = B/A = (Q*R)/A = Q*(R/A) 
%
% Case B is a matrix: X*A = X*(Q*R) = B ==> X = (B/R)*Q' ==> X' = Q*(B/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(A, 2) ~= size(B, 2)) && ~(isa(A, 'double') && isscalar(A)) )
    error('CHEBFUN:CHEBTECH:mrdivide:size', 'Matrix dimensions must agree.');
elseif ( isa(A, 'double') )  % CHEBTECH / double
    if ( ~any(A(:)) )  
        X = B.make(NaN(1, size(B, 2)));
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
elseif ( isa(B, 'double') )  % double / CHEBTECH
    % Here B is a double and A is a CHEBTECH. Do least squares via QR:
    [Q, R] = qr(A, 0);
    
    % Return the transpose for the output.
    X = Q*(B/R).';
elseif ( isa(A, 'chebtech') && isa(B, 'chebtech') )
    error('CHEBFUN:CHEBTECH:mrdivide:funfun', ...
        'Use ./ to divide by a CHEBTECH.');
else
    error('CHEBFUN:CHEBTECH:mrdivide:derp', '%s/%s is not well-defined.', ...
        class(B), class(A));

end

end
