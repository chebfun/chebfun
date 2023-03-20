X = linspace(-1,1,1000); F = tanh(20*X);
subplot(1,2,1)
r = aaa(F,X,'mmax',16); plot(X,F-r(X)), hold on
r = aaa(F,X,'mmax',16,'lawson',50); plot(X,F-r(X)), hold off

Z = exp(1i*pi*linspace(-1,1,1000)); G = exp(Z);
subplot(1,2,2)
r = aaa(G,Z,'tol',1e-6); plot(G-r(Z)), axis equal, hold on
r = aaa(G,Z,'tol',1e-6,'lawson',10); plot(G-r(Z)), axis equal, hold off



