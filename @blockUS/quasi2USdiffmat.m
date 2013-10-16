function L = quasi2USdiffmat(L, dim)

A = blockUS(dim, [-1,1]);
order = size(L,2)-1;

for j = 1:size(L,2)
    %form D^(j-1) term.
    const = feval(L(:,j),0);   % assume constant coeff for now. 

    A = A + const*convert(A,j-1,order)*diff(A, j - 1);
end

end