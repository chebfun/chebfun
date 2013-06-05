function out = poly(f)
%POLY   Polynomial coefficients of a BNDFUN.
%   C = POLY(F) returns the polynomial coefficients of F so that 
%           F(x) = C(1)*x^N + C(2)*x^(N-1) + ... + C(N)*x + C(N+1)
%   
%   Note that unlike the MATLAB POLY command, BNDFUN/POLY can operate on
%   array-valued BNDFUN objects, and hence produce a matrix output. In such
%   instances, the rows of C correspond to the columns of F = [F1, F2, ...].
%   That is, 
%           F1(x) = C(1,1)*x^N + C(1,2)*x^(N-1) + ... + C(1,N)*x + C(1,N+1)
%           F2(x) = C(2,1)*x^N + C(2,2)*x^(N-1) + ... + C(2,N)*x + C(2,N+1).
%   This strange behaviour is a result of MATLAB's decision to return a row
%   vector from the POLY command, even for column vector input.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
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

% Call poly() on the ONEFUN of f
onefunPoly = poly(f.onefun);
n = length(onefunPoly);


% Obtain endpoints of domain, and rescale coefficients if necessary
a = f.domain(1);
b = f.domain(2);

% Convert coefficients on [-1,1] to [a,b].
% TODO: Explain why this formula works/give reference.
if ( a~=-1 || b~=1 )
    % Flip coefficients
    out = onefunPoly(end:-1:1);
    
    % Constants for rescaling
    alpha = 2/(b-a); 
    beta = -(b+a)/(b-a);

    % Rescale coefficients to actual interval
    for j = 0:n-1
        % Need to update coefficients with k>=j.
        k = j:n-1;
        
        % Binomial coefficients, which seem to be more accurate than using
        % MATLAB's NCHOOSEK.
        binom = round(exp(gammaln(k+1)-gammaln(k-j+1)-gammaln(j+1))); 
        
        % Adjust coefficients
        out(j+1) = sum(out(k+1).*binom.*beta.^(k-j).*alpha^(j));
    end
    
    % Flip coefficients back
    out = out(end:-1:1);
end

