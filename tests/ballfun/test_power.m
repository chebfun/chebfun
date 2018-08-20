function pass = test_power( ) 
% Test with function r*cos(lam)*sin(th) to the power 5
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th),[20,21,22]);
g = power(f,5);
h = ballfun(@(r,lam,th)r.^5.*cos(lam).^5.*sin(th).^5,[20,21,22]);

pass(1) = isequal(g,h);
end
