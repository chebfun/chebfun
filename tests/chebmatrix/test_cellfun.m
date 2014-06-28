function pass = test_cellfun(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(666)

% Test that a UniformOutput defaults to false:
r = rand(5,4);
A = chebmatrix(r);
nA = cellfun(@norm, A);
pass(1) = isa(nA, 'chebmatrix') && all(size(nA) == size(A)) && ...
    norm(cell2mat(A.blocks) - abs(r)) == 0;

% Test that a UniformOutput works:
nA = cellfun(@norm, A, 'UniformOutput' ,true);
pass(2) = isnumeric(nA) && all(size(nA) == size(A)) && norm(nA - abs(r)) == 0;

% Test that a function of two variables:
B = chebmatrix(rand(5,4));
AB = cellfun(@times, A, B);
pass(3) = isa(AB, 'chebmatrix') && all(size(AB) == size(A));

% Test that a function of two variables with uniform output:
A = chebmatrix(chebfun(rand(5, 4)));
B = chebmatrix(chebfun(rand(5, 4)));
I = cellfun(@(x,y) sqrt(sum(x.*y)), A, B, 'UniformOutput', true);
pass(4) = isnumeric(I) && all(size(I) == size(A));

end




