function X = mrdivide(A, B)
%/   Right matrix divide for a BNDFUN.
%   A/B divides the BNDFUN A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is a BNDFUN.  Note that in the latter case, formally it
%   is X.' that is returned, as BNDFUN objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case B is a matrix: X*B = A  ==> X = A/B = (Q*R)/B = Q*(R/B)
%
% Case A is a matrix: X*B = X*(Q*R) = A ==> X = (A/R)*Q' ==> X' = Q*(A/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(B, 2) ~= size(A, 2)) && ~(isa(B, 'double') && isscalar(B)) )
    error('CHEBFUN:BNDFUN:mrdivide:size', 'Matrix dimensions must agree.');
    
elseif ( isa(B, 'double') )  % BNDFUN / double
    
    if ( isscalar(B) )
        % Scalar case is easy:
        X = A;                              % Copy A to X
        X.onefun = X.onefun/B;              % mrdivide of the onefun
    else
        % Call MRDIVIDE at the ONEFUN level:
        X = A;
        X.onefun = (A.onefun/B);
        
%         % Alternatively, we could call QR() at the BNDFUN level
%         % For matrix case, we do least squares via QR:
%         [Q, R] = qr(A, 0);
%         X = Q*(R/B);
    end
    
elseif ( isa(A, 'double') )  % double / BNDFUN
    
    % Call MRDIVIDE at the ONEFUN level:
    X = B;
    X.onefun = (A/B.onefun);
    X = X/(.5*diff(B.domain));
    
%     % Alternatively we could call QR() at the BNDFUN level:
%     % Do least squares via QR:
%     [Q, R] = qr(B, 0);
%     % Return the transpose for the output.
%     X = Q*(A/R).';
    
elseif ( isa(B, 'bndfun') && isa(A, 'bndfun') )
    error('CHEBFUN:BNDFUN:mrdivide:bndfunDivBndfun', ...
        'Use ./ to divide BNDFUN by a BNDFUN.');
    
else
    error('CHEBFUN:BNDFUN:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));
    
end

end
