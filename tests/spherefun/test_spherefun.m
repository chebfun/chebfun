function pass = testspherefun( )
% % Main testing file, for now.  
% pass(1) = all(test_constructor( )); 
% pass(2) = all(test_feval( )); 
% pass(3) = all(test_sum2( )); 
% pass(4) = all(test_plus( ));
% pass(5) = all(test_times( )); 
% pass(6) = all(test_power( ));
% pass(7) = all(test_abs( ));
% pass(8) = all(test_diff( ));
% pass(9) = all(test_laplacian( ));
end 



















% function fdf = sph2torus(f,lam,th)
% 
% fdf = real(f(lam,th));
% 
% id = th-pi/2 > 100*eps;
% 
% if ~isempty(id) && any(id(:))
%     fdf(id) = f(lam(id)-pi,pi-th(id));
% end
% 
% end
% 
% function fdf = sphf2cartf(f,lam,th)
% 
% x = cos(lam).*cos(th);
% y = sin(lam).*cos(th);
% z = sin(th);
% 
% fdf = f(x,y,z);
% 
% end


