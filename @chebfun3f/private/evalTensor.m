function T = evalTensor(I, J, K, ff,vectorize)
if vectorize == 1 % we can use the efficient evaluations
n = [numel(I),numel(J),numel(K)];
x = zeros([n(1),1,1]);
x(:,1,1) = I;
X = repmat(x,1,n(2),n(3));
y = zeros([1,n(2),1]);
y(1,:,1) = J;
Y = repmat(y,n(1),1,n(3));
z = zeros([1,1,n(3)]);
z(1,1,:) = K;
Z = repmat(z,n(1),n(2),1);
T = ff(X,Y,Z);
else % we need for loops as f is not vectorizable
T = zeros(size(I,2),size(J,2),size(K,2));
for i = 1:size(I,2)
    for j =1:size(J,2)
        for k = 1:size(K,2)
            T(i,j,k) = ff(I(i),J(j),K(k));
        end
    end
end    
end
end