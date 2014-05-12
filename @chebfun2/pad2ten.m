function [C, V, X, Y] = pad2ten(f, dom)
%CHEBFUN2.PAD2TEN   Evaluate Padua interpolant on a tensor product grid.
%   C = CHEBFUN2.PAD2TEN(F) returns the bivariate Chebyshev coefficients C of
%   the Padua interpolant to the data {X, F}, where X is the Padua grid returned
%   by PADUAPTS(N) for an appropriately chosen value of N.
%
%   [C, V, X, Y] = CHEBFUN2.PAD2TEN(F) returns also the values V of the same
%   interpolant evaluated at an (N+1)x(N+1) point Chebyshev tensor product
%   grid, {X, Y}.
%
%   ... = CHEBFUN2.PAD2TEN(F, [a, b, c, d]) is as above, but when F is given by
%   PADUAPTS(N, [a, b, c, d]).
%
%   Notes: 
%      * The ordering of V is consistent with CHEBFUN2 matrix inputs.
%      * This code is inspired by the algorithm in [1].
%
% See also PADUAPTS.

% References:
%  [1]  Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello
%       "Padua2DM: fast interpolation and cubature at the Padua points in
%       Matlab/Octave."

% Use matrices for small n, DCTs for large n. 
useFFTwhenNisMoreThan = 100; % This decides to switching point.

% Find the appropriate n (so that length(paduapts(n)) == length(f):
m = length(f);
n = round(-1.5 + sqrt(.25+2*m));

% Compute the Padua points:
[x, idx] = paduapts(n);

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
    dct = @(c) flipud(chebtech2.coeffs2vals(flipud(c)));    
    C = dct(dct(G.').');
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

% Evaluate on a tensor product grid:
if ( n < useFFTwhenNisMoreThan )
    % Use Matrices:
    V = Tn1*C*Tn1;   
else
    % Use DCT:
    V = dct(dct(C.').');    
end
% For consistency with Chebfun2.
V = rot90(V, 2); 

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





