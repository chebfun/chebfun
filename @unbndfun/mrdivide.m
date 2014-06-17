function X = mrdivide(A, B)
%/   Right matrix divide for an UNBNDFUN.
%   A/B divides the UNBNDFUN A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is an UNBNDFUN. Note that in the latter case, formally it
%   is X.' that is returned, as UNBNDFUN objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case B is a matrix: X*B = A  ==> X = A/B = (Q*R)/B = Q*(R/B)
%
% Case A is a matrix: X*B = X*(Q*R) = A ==> X = (A/R)*Q' ==> X' = Q*(A/R)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( (size(B, 2) ~= size(A, 2)) && ~(isa(B, 'double') && isscalar(B)) )
    error('CHEBFUN:UNBNDFUN:mrdivide:size', 'Matrix dimensions must agree.');
    
elseif ( isa(B, 'double') )  % UNBNDFUN / double
    
    if ( isscalar(B) )
        % Scalar case is easy:
        X = A;                              % Copy A to X
        X.onefun = X.onefun/B;              % mrdivide of the onefun
    else
        % Call MRDIVIDE at the ONEFUN level:
        X = A;
        X.onefun = (A.onefun/B);
        
%         % Alternatively, we could call QR() at the UNBNDFUN level
%         % For matrix case, we do least squares via QR:
%         [Q, R] = qr(A, 0);
%         X = Q*(R/B);
    end
    
elseif ( isa(A, 'double') )  % double / UNBNDFUN
    error('CHEBFUN:UNBNDFUN:mrdivide:doubleDivUnbndfun', ...
        ['/ does not support the division of a numerical matrix divided by ' ...
        'by an UNBNDFUN.']);
    
elseif ( isa(B, 'unbndfun') && isa(A, 'unbndfun') )
    error('CHEBFUN:UNBNDFUN:mrdivide:unbndfunDivUnbndfun', ...
        'Use ./ to divide UNBNDFUN by an UNBNDFUN.');
    
else
    error('CHEBFUN:UNBNDFUN:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));
    
end

end
