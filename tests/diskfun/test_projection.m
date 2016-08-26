function pass = test_projection( ) 
% Test diskfun projectOntoBMCII() command 

tol = 100*chebfunpref().cheb2Prefs.chebfun2eps;

% Test a constant plus noise
f = diskfun(@(x,y) 1 + 0*x);
% Add noise to the samples then construct a diskfun from the samples
U = sample(f, 16, 17) + (1-2*rand(17,16))*1e-10;
g = diskfun(U);
% Sample the diskfun.
m = 101;
[t,r] = meshgrid(linspace(-pi, pi, m), chebpts(m));
V = g(t, r, 'polar');
% Check that it is constant along the pole

pass(1) = std(V((m-1)/2+1, :)) < tol;
% Check that it has BMC symmetry
id1 = 1:(m-1)/2+1; 
id2 = (m-1)/2+1:m; 
id3 = m:-1:(m-1)/2+1;
pass(2) = norm(V(id1, id1) - V(id3, id2), inf) < tol && ...
          norm(V(id1, id2) - V(id3, id1), inf) < tol;
      
% Test a more complicated function plus noise
f = diskfun(@(x,y) cos(pi*x.*y) + sin(3*pi*x));
% Add noise to the samples then construct a diskfun from the samples
U = sample(f, 16, 17) + (1-2*rand(17,16))*1e-10;
g = diskfun(U);
% Sample the diskfun.
m = 101;
[t,r] = meshgrid(linspace(-pi, pi, m), chebpts(m));
V = g(t,r, 'polar');
% Check that it is constant along the poles
pass(3) = std(V((m-1)/2+1, :)) < tol;
% Check that it has BMC symmetry
id1 = 1:(m-1)/2+1; 
id2 = (m-1)/2+1:m; 
id3 = m:-1:(m-1)/2+1;
pass(4) = norm(V(id1, id1) - V(id3, id2), inf) < 10*tol && ...
          norm(V(id1, id2) - V(id3, id1), inf) < 10*tol;

% Diskfun with a mix of even/pi-period and odd/pi-anti-periodic parts
f = diskfun(@(x,y,z) cos(pi*x.*y) + x.*sin(pi*x.*y));
% Switch the even and odd parts.
temp = f.idxPlus;
f.idxPlus = f.idxMinus;
f.idxMinus = temp;
% The result after projecting f should be exactly zero.
g = projectOntoBMCII(f);
pass(5) = norm(g,inf) < tol;
 
end