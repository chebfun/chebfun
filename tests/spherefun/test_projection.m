function pass = test_projection( ) 
% Test spherefun projectOntoBMCI() command 

tol = 100*chebfunpref().cheb2Prefs.chebfun2eps;

% Test a constant plus noise
f = spherefun(@(x,y,z) 1 + 0*x);
% Add noise to the samples then construct a spherefun from the samples
U = sample(f, 16, 16) + (1-2*rand(16))*1e-10;
g = spherefun(U);
% Sample the spherefun.
m = 101;
[lam,th] = meshgrid(linspace(-pi, pi, m));
V = g(lam, th);
% Check that it is constant along the poles
pass(1) = std(V(end, :)) <= tol;
pass(2) = std(V((m-1)/2+1, :)) < tol;
% Check that it has BMC symmetry
id1 = 1:(m-1)/2+1; 
id2 = (m-1)/2+1:m; 
id3 = m:-1:(m-1)/2+1;
pass(3) = norm(V(id1, id1) - V(id3, id2), inf) < tol && ...
          norm(V(id1, id2) - V(id3, id1), inf) < tol;
      
% Test a more complicated function plus noise
f = spherefun(@(x,y,z) cos(pi*x.*z) + sin(pi*x.*y));
% Add noise to the samples then construct a spherefun from the samples
U = sample(f, 16, 16) + (1-2*rand(16))*1e-10;
g = spherefun(U);
% Sample the spherefun.
m = 101;
[lam,th] = meshgrid(linspace(-pi, pi, m));
V = g(lam, th);
% Check that it is constant along the poles
pass(4) = std(V(end, :)) <= tol;
pass(5) = std(V((m-1)/2+1, :)) < tol;
% Check that it has BMC symmetry
id1 = 1:(m-1)/2+1; 
id2 = (m-1)/2+1:m; 
id3 = m:-1:(m-1)/2+1;
pass(6) = norm(V(id1, id1) - V(id3, id2), inf) < 10*tol && ...
          norm(V(id1, id2) - V(id3, id1), inf) < 10*tol;

% Spherefun with a mix of even/pi-period and odd/pi-anti-periodic parts
f = spherefun(@(x,y,z) cos(pi*x.*z) + x.*sin(pi*x.*y));
% Switch the even and odd parts.
temp = f.idxPlus;
f.idxPlus = f.idxMinus;
f.idxMinus = temp;
% The result after projecting f should be exactly zero.
g = projectOntoBMCI(f);
pass(7) = norm(g,inf) < tol;
 
end