function out = sum(f, dim)
%SUM   Definite integral of a FOURIERTECH on the interval [-PI,PI].
%   SUM(F) is the integral of F from -PI to PI.
%
%   If F is an array-valued FOURIERTECH, then the result is a row vector
%   containing the definite integrals of each column.
%
%   SUM(F, 2) sums over the second dimension of F, i.e., adds up its columns.
%   If F is a scalar-valued FOURIERTECH, this simply returns F.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Get the length of the values:
n = size(f.values, 1);

%%
% Sum across array-valued CHEBTECH columns if dim = 2:
if ( nargin > 1 && dim == 2 )
    f.values = sum(f.values, dim);
    f.coeffs = sum(f.coeffs, dim);
    vscale = max(abs(f.values), [], 1);
    f.epslevel = sum(f.epslevel.*f.vscale, 2)./vscale;
    f.vscale = vscale;
    out = f;
    return
end
   
%%
% Compute the integral if dim = 1:

% Trivial cases:
if ( isempty(f) )    
    % Empty CHEBTECH
    out = []; 
else
    % Trapezium rule:
    out = sum(f.values, 1)*(2*pi/length(f));
end

end
