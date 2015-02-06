function out = deflationFun(Nu, u, r)

% Norm function
normFun = norm(u-r, 'fro');

% Deflator operator
out = Nu/normFun;
end