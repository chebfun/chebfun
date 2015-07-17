function varargout = coeffs2( f ) 
% COEFFS2   Fourier--Chebyshev coefficients of a diskfun 
% 
%  X = COEFFS2( F ) returns the modes of the spherefun in the 
%  Fourier--Chebyshev based. 
% 
% [C, D, R] = COEFFS2( F ) returns a low rank approximation to the modes. 
% 

% Developers note:
% if there was a class between separableApprox and spherefun, then this
% function could implement both coeffs2 for spherefun and diskfun. 

% Calculate the CDR decomposition: 
[C, D, R] = cdr( f ); 

% Find the Fourier coefficients of each slice: 
U = C.coeffs; 
R = R.coeffs;

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