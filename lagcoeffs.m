function c = lagcoeffs(f, n, alp)

if ( nargin < 3 ) 
    alp = 0;
end

d = domain(f);
L = lagpoly(0:n, alp, d);
c = L\f;


end