function pass = test_PTdecomposition( )

% Example 1 :
S = [51,51,51];
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)),S);
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)),S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(1) = isequal(diff(p,[0,1,0]),diff(p2,[0,1,0]))...
       && isequal(diff(p,[0,0,1]),diff(p2,[0,0,1]))...
       && isequal(diff(t,[0,1,0]),diff(t2,[0,1,0]))...
       && isequal(diff(t,[0,0,1]),diff(t2,[0,0,1]));

% Example 2 :
S = [53,53,53];
p = ballfun(@(r,lam,th)cos(r.^2.*sin(th).^2.*cos(lam).*sin(lam)),S);
t = ballfun(@(r,lam,th)sin(r.^2.*sin(th).*cos(th).*sin(lam)),S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(2) = isequal(diff(p,[0,1,0]),diff(p2,[0,1,0]))...
       && isequal(diff(p,[0,0,1]),diff(p2,[0,0,1]))...
       && isequal(diff(t,[0,1,0]),diff(t2,[0,1,0]))...
       && isequal(diff(t,[0,0,1]),diff(t2,[0,0,1]));

% Example 3 :
S = [41,41,41];
p = randnfunball(5,S);
t = randnfunball(5,S);
V = ballfunv.PT2ballfunv(p,t);
[p2, t2] = PTdecomposition(V);
pass(3) = isequal(diff(p,[0,1,0]),diff(p2,[0,1,0]))...
       && isequal(diff(p,[0,0,1]),diff(p2,[0,0,1]))...
       && isequal(diff(t,[0,1,0]),diff(t2,[0,1,0]))...
       && isequal(diff(t,[0,0,1]),diff(t2,[0,0,1]));

if (nargout > 0)
    pass = all(pass(:));
end
end
