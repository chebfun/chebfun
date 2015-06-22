function out = deflationFun(Nu, u, r, p, alp)

% Norm function
normFun = norm(u-r, 'fro')^p;

% Deflator operator
out = Nu*(1/normFun+alp);
end