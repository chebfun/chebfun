function pass = test_constructor( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Can we make a ballfun object:

% Example 1
n = 10; 
f = ballfun( ones(n,n,n) );
g = ballfunv(f,f,f);
pass(1) = 1;

% Example 2
S = [11,12,13];
v = ballfunv(@(x,y,z)x.*z,@(x,y,z)y,@(x,y,z)y.*x,'cart',S);
vx = ballfun(@(x,y,z)x.*z,'cart',S);
vy = ballfun(@(x,y,z)y,'cart',S);
vz = ballfun(@(x,y,z)y.*x,'cart',S);
w = ballfunv(vx,vy,vz);
pass(2) = norm(v-w)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
