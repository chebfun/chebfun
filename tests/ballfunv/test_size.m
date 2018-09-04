function pass = test_size( ) 

% Example 1
Vx = ballfun(ones(20,21,22));
Vy = ballfun(ones(10,3,2));
Vz = ballfun(ones(11,12,1));
S = [20,21,22;10,3,2;11,12,1];
V = ballfunv(Vx,Vy,Vz);
pass(1) = isequal(S,size(V));

if (nargout > 0)
    pass = all(pass(:));
end
end
