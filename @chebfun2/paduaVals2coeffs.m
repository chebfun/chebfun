function [C, V, X, Y] = paduaVals2coeffs(f, dom)
%CHEBFUN2.PADUAVALS2COEFFS   Get Chebyshev coefficients of a Padua interpolant.
%   CHEBFUN2.PADUAVALS2COEFFS(F) returns the bivariate Chebyshev coefficients of
%   the Padua interpolant to the data {X, F}, where X is the Padua grid returned
%   by PADUAPTS(N) for an appropriately chosen value of N.
%
%   [C, V, X, Y] = CHEBFUN2.PADUAVALS2COEFFS(F) returns also the values V of the
%   same interpolant evaluated at an (N+1)x(N+1) point 2nd-kind Chebyshev tensor
%   product grid, {X, Y}.
%
%   ... = CHEBFUN2.PADUAVALS2COEFFS(F, [a, b, c, d]) is as above, but when F is
%   given by PADUAPTS(N, [a, b, c, d]).
%
%   Notes: 
%      * The ordering of C and V is consistent with CHEBFUN2.VALS2COEFFS().
%      * This code is inspired by the algorithm in [1].
%
% See also PADUAPTS, COEFFS2VALS, VALS2COEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% References:
%  [1]  Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello
%       "Padua2DM: fast interpolation and cubature at the Padua points in
%       Matlab/Octave."

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   To compute the values on the Padua interpolant on a tensor product grid, we
%   _could_ use V = COEFFS2VALS(PADUAVALS2COEFFS(F)), however it's slightly more
%   efficient to use [~, V] = PADUAVALS2COEFFS(F) as the DCT matrices required
%   for COEFFS2VALS() have already been constructed in PADUAVALS2COEFFS().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use matrices for small n, DCTs for large n. 
useFFTwhenNisMoreThan = 100; % This decides to switching point.

% Find the appropriate n (so that length(paduapts(n)) == length(f):
m = length(f);
n = round(-1.5 + sqrt(.25+2*m));

% Compute the Padua points:
[x, idx] = paduapts(n);

%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHEBYSHEV COEFFICIENTS  %%%%%%%%%%%%%%%%%%%%%%%%

% Interpolation weights:
w = 0*x(:,1) + 1./(n*(n+1));
idx1 = all(abs(x) == 1, 2);
w(idx1) = .5*w(idx1);
idx2 = all(abs(x) ~= 1, 2);
w(idx2) = 2*w(idx2);

% Make G: (idx is returned by PADUAPTS()).
G = zeros(size(idx));
G(idx) = 4*w.*f;

% Compute bivariate coefficients:
if ( n < useFFTwhenNisMoreThan )
    % Use Matrices:
    Tn1 = cos((0:n).'*(0:n)*pi/n);
    Tn2 = cos((0:n+1).'*(0:n+1)*pi/(n+1));
    C = Tn2*G*Tn1;
else
    % Use DCT:
    dct = @(c) chebtech2.coeffs2vals(c);    
    C = rot90(dct(dct(G.').'), 2);
end
% Modify a few entries:
C(1,:) = .5*C(1,:);
C(:,1) = .5*C(:,1);
C(1,end) = .5*C(1,end);
C(end,:) = [];

% Take upper-left triangular part: (equivalent to C = fliplr(triu(fliplr(C))))
C = triu(C(:,end:-1:1));
C = C(:,end:-1:1);  

if ( nargout < 2 )
    % No need to go any further!
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATE ON A TENSOR GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate on a tensor product grid:
if ( n < useFFTwhenNisMoreThan )
    % Use Matrices:
    V = Tn1*C*Tn1;   
else
    % Use DCT:
    V = dct(dct(C.').');    
end

if ( nargout < 3 )
    % No need to go any further!
    return
end

% Return the tensor product grid:
if ( nargin < 2 )
    dom = [-1, 1, -1, 1];
end
xnp1 = chebpts(n+1, dom(1:2));
ynp1 = chebpts(n+1, dom(3:4));
[X, Y] = meshgrid(xnp1, ynp1);

end
