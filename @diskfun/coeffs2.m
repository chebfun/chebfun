function varargout = coeffs2( f, m, n ) 
% COEFFS2   Fourier--Chebyshev coefficients of a diskfun 
% 
%  X = COEFFS2( F ) returns the modes of the diskfun in the 
%  Fourier--Chebyshev bases. 
% 
% [C, D, R] = COEFFS2( F ) returns a low rank approximation to the modes. 
% 

% Developers note:
% if there was a class between separableApprox and diskfun, then this
% function could implement both coeffs2 for diskfun and diskfun. 

% Calculate the CDR decomposition: 
[C, D, R] = cdr( f ); 

if nargin == 1
    % Find the  coefficients of each slice: 
    U = C.coeffs; 
    R = R.coeffs;
else
    if nargin == 2
        n = m;
    end
    % Find the  coefficients of each slice: 
    U = chebtech2.alias(C.coeffs,n); 
    R = trigtech.alias(R.coeffs,m);
end

% Prepare the output. Keep in low rank form if nargin > 1.
if ( nargout <= 1 ) 
    varargout = {U*D*R.'};
elseif ( nargout <= 3 )
    varargout = {U, D, R};
else
    error('DISKFUN::COEFFS::NARGOUT',...
            'Too many output arguments')
end

end