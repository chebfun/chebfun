function out = poly(f)
%POLY   Polynomial coefficients of a CHEBTECH.
%   C = POLY(F) returns the polynomial coefficients of F so that 
%           F(x) = C(1)*x^N + C(2)*x^(N-1) + ... + C(N)*x + C(N+1)
%   
%   Note that unlike the MATLAB POLY command, CHEBTECH/POLY can operate on
%   array-valued CHEBTECH objects, and hence produce a matrix output. In such
%   instances, the rows of C correspond to the columns of F = [F1, F2, ...].
%   That is, 
%           F1(x) = C(1,1)*x^N + C(1,2)*x^(N-1) + ... + C(1,N)*x + C(1,N+1)
%           F2(x) = C(2,1)*x^N + C(2,2)*x^(N-1) + ... + C(2,N)*x + C(2,N+1).
%   This strange behaviour is a result of MATLAB's decision to return a row
%   vector from the POLY command, even for column vector input.
%
% See also CHEBCOEFFS, FOURCOEFFS, LEGCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mathematical reference]: Section 3.3 Mason & Handscomb, "Chebyshev
% Polynomials". Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

% Flip the Chebyshev coefficients to match Matlab ordering:
coeffs = flipud(f.coeffs);
[n, m] = size(coeffs);

% Coefficients on the unit interval:
if ( n == 1 )
    % Constant case:
    out = coeffs;
    
elseif ( n == 2 )
    % Linear case:
    out = coeffs.';
    
else
    % General case:
    
    % Initialise working vectors:
    tn = zeros(m, n);
    tnold1 = tn;
    tnold1(:,2) = 1;
    tnold2 = tn;
    tnold2(:,1) = 1;
    zer = zeros(m, 1);
    
    % Transpose coeffs (as output is in this form):
    coeffs = coeffs.';
    
    % Initialise output:
    out = zeros(m, n);
    
    % Initial step:
    out(:,[1,2]) = [ zer coeffs(:,end).*tnold2(:,1) ] + ...
        bsxfun(@times, coeffs(:,end-1), tnold1(:,[2,1]));

    % Recurrence:
    for k = 3:n
        tn(:,1:k)  = [ zer 2*tnold1(:,1:k-1) ] - [ tnold2(:,1:k-2), zer, zer ];
        out(:,1:k) = bsxfun(@times, coeffs(:,end-k+1), tn(:,k:-1:1)) + ...
                        [ zer, out(:,1:k-1) ];                 
        tnold2 = tnold1;
        tnold1 = tn;
    end
    
end

end