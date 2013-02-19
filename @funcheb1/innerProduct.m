function out = innerProduct(f, g)
%INNERPRODUCT
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   FUNCHEB1 objects F and G.

if ( ~isa(f, 'funcheb1') || ~isa(g, 'funcheb1') )
    error('CHEBFUN:FUNCHEB1:InnerProduct:input', ...
        'innerProduct() only operates on two FUNCHEB1 objects.');
end

% Prolong to sum of current lengths (so that quadrature is exact):
n = length(f) + length(g);
f = prolong(f, n);
g = prolong(g, n);

% Compute Clenshaw-Curtis quadrature weights:
w = funcheb1.quadwts(n);

% Compute the inner product via a weighted discrete inner product:
out = (spdiags(w.', 0, n, n) * f.values)' * g.values;

end