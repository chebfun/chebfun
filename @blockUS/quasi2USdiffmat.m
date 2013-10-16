function L = quasi2USdiffmat(L, dim)

A = blockUS(dim, [-1, 1]);


if ( isa(L, 'blockCoeff') )
    c = L.coeffs;
else
    c = L;
end
c = fliplr(c);
order = numel(c) - 1;

L = 0*speye(dim);
for j = 1:size(c, 2)
    %form D^(j-1) term.
    const = feval(c{j}, 0);   % assume constant coeff for now. 
    L = L + const*convert(A, j, order)*diff(A, j - 1);
end

end