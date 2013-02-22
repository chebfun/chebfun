function out = innerProduct(f, g)
%INNERPRODUCT  Compute the inner product of two FUNCHEB objects.
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   FUNCHEB objects F and G (conjugate linear in F).
%
%   If F and/or G are vector-valued FUNCHEB objects, then the result is a matrix
%   whose i-j entry is the inner product of the ith column of F with the jth
%   column of G.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'funcheb') || ~isa(g, 'funcheb') )
    error('CHEBFUN:FUNCHEB:InnerProduct:input', ...
        'innerProduct() only operates on two FUNCHEB objects.');
end

% Prolong to sum of current lengths (so that quadrature is exact):
n = length(f) + length(g);
f = prolong(f, n);
g = prolong(g, n);

% Compute Clenshaw-Curtis quadrature weights:
w = f.quadwts(n);

% Compute the inner product via a weighted discrete inner product:
out = bsxfun(@times, w.', f.values)' * g.values;

end