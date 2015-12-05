function sample_error = test_sampleError( h, g ) 
m = 6; n = m;  
[x, y] = getPoints( m, n ); 
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
% hn = pi/length(g.cols);  % Have to adjust for the shift in y points. Ugly!
% approx = feval(g.cols,(y-hn)/pi-.5) * g.blockDiag * feval(g.rows,x/pi)';
% sample_error = norm( F - approx , inf );
% [C,D,R] = cdr(g);
% approx = feval(g.cols,y/pi) * g.blockDiag * feval(g.rows,x/pi)';
approx = fevalm(g,x,y);
sample_error = norm( F(:) - approx(:) , inf );
end