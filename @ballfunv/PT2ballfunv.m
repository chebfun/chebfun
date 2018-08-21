function v = PT2ballfunv(P,T)
% Take 2 scalar fields P and T  in Cheb-Fourier-Fourier basis and return the vector field
% v = curl(curl(rP)) + curl(rT)
%
% Also see PTDECOMPOSITION

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

v = curl(ballfunv.rcurl(P)) + ballfunv.rcurl(T);
end