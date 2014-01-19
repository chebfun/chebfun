function pass = chebfun2v_subsref
% Test for Chebfun2v subsref. 
% Alex Townsend, March 2013. 

pass = 1;
f = @(x,y) cos(x); g=chebfun2v(f,f);  % any fun2.
c = chebfun(f); d = c; 
try
    % subrefs working with single reference
    fsub = g.xcheb;
    % get working
    fget = get(g,'xcheb');
    % hard code subreferencing. Take one dimension slices. 
%     pass(1) = (norm(g(:,pi/6)-c)<1e-14);
%     pass(2) = (norm(g(pi/6,:)-f(pi/6))<1e-14);

    pass = all(pass); 
catch
    pass = 0;
end
end