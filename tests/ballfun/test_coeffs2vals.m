function pass = test_coeffs2vals( ) 
% Test with function sin(r)*th*lam
f = ballfun(@(r,lam,th)sin(r).*lam.*th,[20,20,20]);

pass(1) = (max(max(max(abs(ballfun.vals2coeffs(ballfun.coeffs2vals(f.coeffs))-f.coeffs))))<1e-15);

end
