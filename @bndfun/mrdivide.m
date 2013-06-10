function X = mrdivide(A, B)
%/   Right matrix divide for a BNDFUN.
%
%   A/B divides the BNDFUN A by a scalar B. More generally, it gives the
%   least-squares solution (with respect to the continuous L^2 norm) to X*B = A
%   when either A or B is a BNDFUN.  Note that in the latter case, formally
%   it is X.' that is returned, as BNDFUN objects are always columns.
%
% See also QR, RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
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
        % For matrix case, we do least squares via QR:
        [Q, R] = qr(A, 0);
        X = Q*(R/B);
    end
elseif ( isa(A, 'double') )  % double / BNDFUN
    % Do least squares via QR:
    [Q, R] = qr(B, 0);
    
    % Return the transpose for the output.
    X = Q*(A/R).';
    
    % TODO: Do we want to compute MRDIVIDE using the method above, or call
    % MRDIVIDE  at the ONEFUN level, e.g. via
    %     X = B;
    %     X.onefun = A/B.onefun*.5*diff(B.domain);
    
elseif ( isa(B, 'BNDFUN') && isa(A, 'BNDFUN') )
    error('CHEBFUN:BNDFUN:mrdivide:BndfunDivBndfun', ...
        'Use ./ to divide by a BNDFUN.');
else
    error('CHEBFUN:BNDFUN:mrdivide:badArg', '%s/%s is not well-defined.', ...
        class(A), class(B));
    
end

end
