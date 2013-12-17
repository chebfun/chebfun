function f = restrict(f, dom)
% RESTRICT  Restrict the domain of a chebfun2.

f.cols = restrict(f.cols, dom(3:4)); 
f.rows = restrict(f.rows, dom(1:2));
f.domain = dom;

%f = simplify( f ); 
end