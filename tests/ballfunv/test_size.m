function pass = test_size( ) 

% Example 1
Vx = ballfun(@(x,y,z)1);
Vy = ballfun(@(x,y,z)1);
Vz = ballfun(@(x,y,z)x);
S = [1,1,1;1,1,1;2,3,3];
V = ballfunv(Vx,Vy,Vz);
pass(1) = isequal(S,size(V));

if (nargout > 0)
    pass = all(pass(:));
end
end
