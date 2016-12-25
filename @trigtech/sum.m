function out = sum(f, dim)
%SUM   Definite integral of a TRIGTECH on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.
%
%   If F is an array-valued TRIGTECH, then the result is a row vector
%   containing the definite integrals of each column.
%
%   SUM(F, 2) sums over the second dimension of F, i.e., adds up its columns.
%   If F is a scalar-valued TRIGTECH, this simply returns F.
%
% See also CUMSUM, DIFF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the length of the values:
n = size(f.values, 1);

%%
% Sum across array-valued TRIGTECH columns if dim = 2:
if ( nargin > 1 && dim == 2 )
    f.values = sum(f.values, dim);
    f.coeffs = sum(f.coeffs, dim);
    f.isReal = all(f.isReal);
    out = f;
    return
end
   
%%
% Compute the integral if dim = 1:

% Trivial cases:
if ( isempty(f) )    
    % Empty TRIGTECH
    out = []; 
else
    % By orthogonality of complex exponentials integral is 2x 
    % zeroth coefficient: 
    out = 2*f.coeffs( floor((n+2)/2), :); 
end

% Return a real result if f is real:
out(:,f.isReal) = real(out(:,f.isReal));

end
