function y = feval(f, x)
%FEVAL   Evaluate a TRIGTECH.
%   Y = FEVAL(F, X) Evaluation of the TRIGTECH F at points X via
%   Horner's scheme.
%
%   If size(F, 2) > 1 then FEVAL returns values in the form [F_1(X), F_2(X),
%   ...], where size(F_k(X)) = size(X).
%
%   Example:
%     f = trigtech(@(x) exp(cos(pi*x)) );
%     x = linspace(-1, 1, 1000);
%     fx = feval(f, x);
%     plot(x,fx,'r-')
%

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    y = [];
    return 
end

% Reshape x to be a column vector for passing to polyval:
[ignore, m] = size(f);
% c = f.coeffs(end:-1:1, :);
% N = size(fourierCoeff, 1);
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

% Evaluate using some variant of Horner's scheme
y = f.horner(x, f.coeffs, all(f.isReal));

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex(1), m*numel(x)/sizex(1));
end

end
