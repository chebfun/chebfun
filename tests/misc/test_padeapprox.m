% Test for padeapprox.m.

function pass = test_padeapprox(pref)

tol = 1e-10;

f = @(x) (x.^4 - 3)./((x + 3.2).*(x - 2.2));
[r, a, b, mu, nu, poles, residues] = padeapprox(f, 10, 10);
pass(1) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(poles) - [-3.2 ; 2.2])) < tol);

c = [1 1i];
[r, a, b] = padeapprox(c, 0, 1);
pass(2) = (abs(a-1) < 1e-10) && (abs(b(2)-(-1i)) < tol); 

f = @(x) x./(1-x);
[r, a, b, mu, nu, poles, residues] = padeapprox(f, 5, 6, [], 0.5);
pass(3) = (mu == 1) && (nu == 1) && norm(a-[0;1])<tol && ...
    norm(b-[1;-1])<tol && abs(poles-1)<tol && abs(residues+1)<tol;


end
