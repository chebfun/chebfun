function pass = test_mean( ) 
% Test with function 1
f = ballfun(@(r,lam,th)1,[20,20,20]);
I = mean(f);
pass(1) = (abs(I-1)<1e-15);
end
