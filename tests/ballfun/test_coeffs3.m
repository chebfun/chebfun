function pass = test_coeffs3( ) 

% Example 1
f = ballfun(@(r,lam,th)1);
F = coeffs3(f,3,4,5);
exact = zeros(3,4,5);
exact(1,3,3) = 1;
pass(1) = isequal(F,exact);

% Example 2
exact = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
F = coeffs3(exact,50,51,52);
f = ballfun(F,"coeffs");
pass(2) = isequal(f,exact) && isequal(size(F),[50,51,52]);

% Example 3
exact = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
F = coeffs3(exact);
f = ballfun(F,"coeffs");
pass(3) = isequal(f,exact);

% Example 4
exact = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
F = coeffs3(exact,50);
f = ballfun(F,"coeffs");
pass(4) = isequal(f,exact) && isequal(size(F),[50,50,50]);

% Example 5
exact = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
F = coeffs3(exact,50,51);
f = ballfun(F,"coeffs");
pass(5) = isequal(f,exact) && isequal(size(F),[50,51,51]);

if (nargout > 0)
    pass = all(pass(:));
end
end