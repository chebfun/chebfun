function varargout = coeffs2( f ) 
% COEFFS2   Fourier coefficients of a spherefun 
% 
%  X = COEFFS2( F ) returns the 2D Fourier modes of the spherefun, viewed
%  as a doubly periodic function. 
% 
% [C, D, R] = COEFFS2( F ) returns a low rank approximation to the 2D
% Fourier modes. 
% 

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
    error('SPHEREFUN::COEFFS::NARGOUT',...
            'Too many output arguments')
end

end