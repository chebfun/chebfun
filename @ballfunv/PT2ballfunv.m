function v = PT2ballfunv(P,T)
% Take 2 scalar fields P and T  in Cheb-Fourier-Fourier basis and return the vector field
% v = curl(curl(rP)) + curl(rT)
v = curl(ballfunv.rcurl(P)) + ballfunv.rcurl(T);
end