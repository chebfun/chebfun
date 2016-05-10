function cU = chebTcoeffs2chebUcoeffs(cT)
%CHEBTCOEFFS2CHEBUCOEFFS   Convert Chebyshev-T coeffs. to Chebyshev-U coeffs.
%   CU = CHEBTCOEFFS2CHEBUCOEFFS(CT) takes a column vector of coefficients of a
%   polynomial represented as an expansion in the Chebyshev polynomials of the
%   first kind T_k(x) (ordered starting with the coefficient of T_0(x)) and
%   returns the column vector of coefficients of the same polynomial expanded
%   in Chebyshev polynomials of the second kind U_k(x) (ordered starting with
%   the coefficient of U_0(x)).
%
%   If CT is a matrix, the conversion is performed column-wise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(cT) )
    cU = [];
    return;
end

% We compute 2nd-kind coefficients from the 1st-kind coefficients
% using the recurrence T_n(x) = (1/2)*(U_n(x) - U_{n-2}(x)).  The
% coefficient of U_n requires the coefficients of T_n and T_{n + 2}, so we
% pad with two zero-coefficients initially.
nCols = size(cT, 2);
cU = [cT ; zeros(2,nCols)];

% Run the recurrence.
cU(1,:) = 2*cT(1,:);
cU = 0.5*[cU(1:end-2,:) - cU(3:end,:) ; cU(end-1:end,:)];
cU = cU(1:end-2,:);

end
