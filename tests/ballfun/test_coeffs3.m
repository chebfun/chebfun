function pass = test_coeffs3( pref ) 

% Example 1
f = ballfun(@(r,lam,th)1);
F = coeffs3(f,3,4,5);
exact = zeros(3,4,5);
exact(1,3,3) = 1;
pass(1) = isequal(F,exact);

% Example 2
exact = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2));
F = coeffs3(exact,50,50,50);
f = ballfun(F,"coeffs");
pass(2) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end