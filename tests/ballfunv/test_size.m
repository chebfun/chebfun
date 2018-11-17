function pass = test_size( ) 

% Example 1
Vx = ballfun(ones(21,20,22));
Vy = ballfun(ones(11,2,2));
Vz = ballfun(ones(11,12,4));
S = [1,1,1;1,1,1;1,1,1];
V = ballfunv(Vx,Vy,Vz);
pass(1) = isequal(S,size(V));

if (nargout > 0)
    pass = all(pass(:));
end
end
