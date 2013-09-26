function S = sum(A)

n = dim(A);
d = domain(A);
numint = length(d)-1;
S = cell(1,numint);
for k = 1:numint
    S{k} = chebtech2.quadwts(n(k));
end
S = cell2num(S);
S = S * (d(end)-d(1))/2;

end
