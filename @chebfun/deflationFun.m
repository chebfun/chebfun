function out = deflationFun(Nu, u, r, p, alp)

if (isa(r,'chebfun'))
    r=  mat2cell(r);
end
% Norm function
normFun = 1;
for rCounter = 1:length(r)
    normFun = normFun*norm(u-r{rCounter}, 'fro')^p;
end

% Deflator operator
out = Nu*(1/normFun+alp);
end