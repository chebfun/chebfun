function L = quasi2USdiffmat(L, dim)

A = blockUS(dim, L.domain);

if ( isa(L, 'blockCoeff') )
    % TODO: Still possible?
    c = L.coeffs;
else
    c = L.coeffForm;
end
c = fliplr(c);
order = numel(c) - 1;

L = 0*speye(sum(dim));
for j = 1:size(c, 2)
    %form D^(j-1) term.
    %const = feval(c{j}, 0);   % assume constant coeff for now. 
    L = L + convert(A, j-1, order-1)*mult(A, c{j}, j-1)*diff(A, j - 1);
end

% if ( dim < 200 ) 
%     L = full(L); 
% end


end