function N = mtimes(N,u)
%MTIMES 

if isa(u,'chebfun2')
    op = N.op;
    N =op(u);     
elseif isa(N,'chebop2') && isa(u,'double')
    N.coeffs = u*N.coeffs; 
    op = N.op; 
    N.op = @(v) u*op(v);    
elseif isa(u,'chebop2') && isa(N,'double')
    N = mtimes(u,N);
else
    error('CHEBOP2:MTIMES','Can only times a chebop2 by a double or forward apply to a chebfun2.'); 
end

end