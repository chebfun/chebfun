function pass = test_uplus( ) 
% Test with function rand = +rand
S = [19,12,27];
F1 = ballfun(@(r,lam,th)r.*cos(lam),S);
F2 = ballfun(@(r,lam,th)r.*sin(lam).*sin(th),S);
F3 = ballfun(@(r,lam,th)cos(lam).*cos(th),S);
F = ballfunv(F1,F2,F3);
G = +F;
pass(1) = isequal(F,G);
end
