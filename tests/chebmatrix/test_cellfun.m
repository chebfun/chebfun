function pass = test_cellfun(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(666)

%% Test that a UniformOutput defaults to false:
R1 = rand(5,4);
A = chebmatrix(R1);
nA = cellfun(@norm, A);
pass(1) = isa(nA, 'chebmatrix') && all(size(nA) == size(A)) && ...
    ( norm(cell2mat(A.blocks) - abs(R1)) == 0 );

%% Test that a UniformOutput works:
nA = cellfun(@norm, A, 'UniformOutput' ,true);
pass(2) = isnumeric(nA) && all(size(nA) == size(A)) && (norm(nA - abs(R1)) == 0);

%% Test that a function of two variables:
R2 = rand(5,4);
B = chebmatrix(R2);
AB = cellfun(@times, A, B);
pass(3) = isa(AB, 'chebmatrix') && all(size(AB) == size(A))&& ...
    (norm(cell2mat(AB.blocks) - R1.*R2) == 0);

%% Test that a function of two variables with uniform output:
A = chebmatrix(chebfun(R1));
B = chebmatrix(chebfun(R2));
I = cellfun(@(x,y) sqrt(sum(x.*y)), A, B, 'UniformOutput', true);
% Compute the correct answer in a non-slick way to compare against:
corr = zeros(1,4);
for colCounter = 1:size(A, 2)
    corr(colCounter) = sqrt(sum(A{colCounter}.*B{colCounter}));
end

pass(4) = isnumeric(I) && all(size(I) == size(A)) && (norm(I-corr) == 0);

end




