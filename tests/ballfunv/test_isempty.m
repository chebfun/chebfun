function pass = test_isempty( ) 

% Example 1:
v = ballfunv();
pass(1) = isempty(v);

% Example 2:
f = ballfun(1);
g = ballfun();
v = ballfunv(f,g,f);
pass(2) = ~isempty(v);

% Example 3:
f = ballfun(1);
v = ballfunv(f,f,f);
pass(3) = ~isempty(v);

if (nargout > 0)
    pass = all(pass(:));
end
end
