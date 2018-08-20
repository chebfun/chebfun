function pass = test_plus( ) 
S = [19,20,21];
% Example 1
f = ballfun(ones(20,21,22));
V1 = ballfun.coeffs2vals(f.coeffs);
g = f+f;
V2 = ballfun.coeffs2vals(g.coeffs);
pass(1) = (max(max(max(abs(V2-2*V1))))==0);

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(th)+1,S);
g = f+1;
exact = ballfun(@(r,lam,th)r.^2.*cos(th)+2,S);
pass(2) = isequal(g,exact);

% Example 3
f = ballfun(@(r,lam,th)r.*sin(th)+1,S);
g = 3+f;
exact = ballfun(@(r,lam,th)r.*sin(th)+4,S);
pass(3) = isequal(g,exact);

% Example 4
f = ballfun(@(r,lam,th)r.*sin(th).^2.*cos(lam)+1,S);
g = 3+f+2;
exact = ballfun(@(r,lam,th)r.*sin(th).^2.*cos(lam)+6,S);
pass(4) = isequal(g,exact);

if (nargout > 0)
    pass = all(pass(:));
end

end
