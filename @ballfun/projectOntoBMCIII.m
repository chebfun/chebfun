function f = projectOntoBMCIII( f )
% PROJECTONTOBMCI  Projection onto BMC-III symmetry.
%
% f = projectOntoBMCIII(f) is the orthogonal projection of f onto BMC-III
% symmetry, i.e., a function that is
% 1. even in theta for every even wave number in lambda;
% 2. odd in theta for every odd wave number in lambda;
% Additionally, for all but k=0 wavenumber lambda the resulting projection
% enforces the ballfun is zero at the poles. 
%
% The projection is orthogonal, i.e., the correction matrix to fix up the
% structure has the smallest possible Frobenius norm.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
   return
end

% Get the tensor of coefficients
F = f.coeffs;

% Permute F
F = permute(F,[1,3,2]);

% Loop over lambda
for l = 1:size(F,3)
    % Update the matrix
    F(:,:,l) = projectOntoRTheta(F(:,:,l));
end

% Permute back
F = permute(F,[1,3,2]);
f = ballfun(F, 'coeffs');
end

function X = projectOntoRTheta( X )
% Project the matrix of Chebyshev--Fourier coefficients onto a BMC-II
% function. The projection is orthogonal, i.e., the correction matrix to
% fix up the structure has the smallest possible Frobenius norm. 

% Get the discretization
[m,n] = size(X);

zeroMode = floor(n/2)+1;
evenModes = [fliplr(zeroMode-2:-2:1) zeroMode+2:2:n];  % Not including the zero mode
oddModes = [fliplr(zeroMode-1:-2:1) zeroMode+1:2:n];

% First do the zero-periodic mode in theta. 
% Enforce the expansion is even in r.
I = eye( m ); A = I(2:2:end,:); 
C = A \ ( A * X(:,zeroMode) ); 
C = I \ C; 

% Update coeff matrix: 
X(:,zeroMode) = X(:,zeroMode) - C; 

% Second do the even-periodic, non-zero modes in theta.
% Enforce these are zero at the pole and that the expansion is even in 
% theta
% Vectors [ T_k(0) ]: 
A = real( (-1i).^(0:m-1) );   % A * X should be the zero vector. 

% Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
% norm: 
C = A \ ( A * X(:,evenModes) ); 
C = I \ C; 
X(:,evenModes) = X(:,evenModes) - C;

% Now project onto odd-anti-periodic.  Nothing special has to be done at
% the poles since enforcing the expansion is odd in r guarantees that it
% sums to zero.
A = I(1:2:end,:); 
C = A \ ( A * X(:,oddModes) ); 
C = I \ C; 
X(:,oddModes) = X(:,oddModes) - C; 
end