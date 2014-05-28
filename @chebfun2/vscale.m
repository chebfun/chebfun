function vscl = vscale( f ) 

[m, n] = length( f ); 
m = max( m, 9 ); 
n = max( n, 9 ); 
vals = chebpolyval2( f, m, n); 
vscl = max( abs( vals(:) ) ); 

end