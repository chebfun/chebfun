function op = singOp2SmoothOp( op, exponents, tol )
%convert a singular opreator to a smooth operator.
if ( all(abs(exponents) > 100*tol ) )
    % both exponents are non trivial
    op = @(x) op(x)./((1+x).^(exponents(1)).*(1-x).^(exponents(2)));
elseif ( abs(exponents(1)) > 100*tol )
    % left only
    op = @(x) op(x)./(1+x).^(exponents(1));
elseif ( abs(exponents(2)) > 100*tol )
    % right only
    op = @(x) op(x)./(1-x).^(exponents(2));
end