% Test file for EXPINTEG/PHIFUN:

function pass = test_phiFun()

tol = 1e-12;

% First four phi functions (exact):
phi0ex = @(z) exp(z);
phi1ex = @(z) (exp(z) - 1)./z;
phi2ex = @(z) (exp(z) - z - 1)./z.^2;
phi3ex = @(z) (exp(z) - z.^2/2 - z - 1)./z.^3;

% Get the same functions with PHIFUN:
phi0 = expinteg.phiFun(0);
phi1 = expinteg.phiFun(1);
phi2 = expinteg.phiFun(2);
phi3 = expinteg.phiFun(3);

% Grid for comparisons:
[xx, yy] = meshgrid([-10:-1, 1:10]);                           
zz = xx + 1i*yy;

% Compare:
pass(1) = max(max(abs(phi0(zz) - phi0ex(zz)))) < tol;
pass(2) = max(max(abs(phi1(zz) - phi1ex(zz)))) < tol;
pass(3) = max(max(abs(phi2(zz) - phi2ex(zz)))) < tol;
pass(4) = max(max(abs(phi3(zz) - phi3ex(zz)))) < tol;

end
