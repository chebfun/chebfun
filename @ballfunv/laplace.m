function f = laplace(v)
% LAPLACE Laplacian of a BALLFUNV
%   LAPLACE(v) is the laplacian of the BALLFUNV v
f = grad(div(v))-curl(curl(v));
end
