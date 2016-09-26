function out = innerProduct(f, g)
%INNERPRODUCT   Compute the inner product of two TRIGTECH objects.
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   TRIGTECH objects F and G (conjugate linear in F).
%
%   If F and/or G are array-valued TRIGTECH objects, then the result is a matrix
%   whose i,j entry is the inner product of the ith column of F with the jth
%   column of G.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'trigtech') || ~isa(g, 'trigtech') )
    error('CHEBFUN:TRIGTECH:innerProduct:input', ...
        'innerProduct() only operates on two TRIGTECH objects.');
end

% Prolong to sum of current lengths (so that quadrature is exact):
n = length(f) + length(g);
f = prolong(f, n);
g = prolong(g, n);

% Compute trapezium quadrature weights:
w = f.quadwts(n);

% Compute the inner product via a weighted discrete inner product:
out = bsxfun(@times, w.', f.values)' * g.values;

% Force real output if the inputs are real:
id = logical(double(f.isReal)' * double(g.isReal));
out(id) = real(out(id));

% Force non-negative output if the inputs are equal:
if ( isequal(f, g) )
    dout = diag(diag(out));
    out = out - dout + abs(dout);
end

end
