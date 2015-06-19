function out = deflationFun(Nu, u, r, alp)

% Norm function
normFun = norm(u-r, 'fro')^2;

% Deflator operator
out = Nu*(1/normFun+alp);
end