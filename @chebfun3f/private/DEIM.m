function indices = DEIM(U)
indices = [];
[~, I] = max(abs(U(:,1)));
indices = [indices,I];
for l = 2:size(U,2)
    c = U(indices,1:(l-1)) \ U(indices,l);
    r = U(:,l) - U(:,1:(l-1))*c;
    [~, I] = max(abs(r));
    indices = [indices,I];
end
end