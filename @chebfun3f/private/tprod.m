function X = tprod(X,U,V,W)
% computes X x_1 U x_2 V x_3 W
n = size(X);
m = [size(U,1),size(V,1),size(W,1)];
X = reshape(U*reshape(X,[n(1),n(2)*n(3)]),[m(1),n(2),n(3)]);
X = permute(reshape(V*reshape(permute(X,[2,1,3]),[n(2),m(1)*n(3)]),[m(2),m(1),n(3)]),[2,1,3]);
X = permute(reshape(W*reshape(permute(X,[3,2,1]),[n(3),m(2)*m(1)]),[m(3),m(2),m(1)]),[3,2,1]);
end
