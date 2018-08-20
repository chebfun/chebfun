function pass = test_uplus( ) 
% Test with function rand = +rand
S = [22,23,14];
f = cheb.galleryballfun('random',S);
g = +f;

pass(1) = isequal(f,g);
end
