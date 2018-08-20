function pass = test_uminus( ) 
% Test with function rand : rand + -rand = 0
S = [20,21,22];
f = cheb.galleryballfun('random',S);
g = -f;

pass(1) = iszero(f+g);
end
