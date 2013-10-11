function S = sum(A)

n = dim(A);
d = A.domain;

numint = length(d)-1;
S = cell(1,numint);
for k = 1:numint
    S{k} = chebtech2.quadwts(n(k)) * (d(k+1)-d(k))/2;
end
S = cat(2,S{:});

end
