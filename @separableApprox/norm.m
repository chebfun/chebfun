function [normF, normloc] = norm( f, p )
%NORM       Norm of a SEPARABLEAPPROX object.
% For SEPARABLEAPPROX objects:
%    NORM(F) = sqrt(integral of abs(F)^2).
%    NORM(F, 2) = largest singular value of F.
%    NORM(F,'fro') is the same as NORM(F).
%    NORM(F,'nuc') = sum of singular values of F.
%    NORM(F, 1) = NOT IMPLEMENTED.
%    NORM(F, inf) = global maximum in absolute value.
%    NORM(F, 'max') = global maximum in absolute value.
%    NORM(F, 'min') = NOT IMPLEMENTED.
%
% Furthermore, the inf norm for SEPARABLEAPPROX objects also returns a 
% second output, giving a position where the max occurs.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to Frobenius norm.
    p = 'fro';
end

if ( isempty( f ) )  
    % Empty chebfun has norm 0.
    normF = [];
    
else
    switch ( p )  % Different cases on different norms.
        case 1
            error('CHEBFUN:SEPARABLEAPPROX:norm:norm', ...
                'SEPARABLEAPPROX does not support L1-norm, yet.');

        case {'fro'}  % Definite integral of f.^2
            % L^2-norm is sum of squares of sv.
            normF = sqrt( sum( svd( f ).^2 ) );  
            
        case {inf, 'inf', 'max'}            
            if ( isreal(f) )
                [Y, X] = minandmax2(f);
                [normF, idx] = max(abs(Y));
                normloc = X(idx, :);
            else
                [Y, X] = minandmax2(conj(f).*f);
                [normF, idx] = max(sqrt(abs(Y)));
                normloc = X(idx, :);
            end
            
        case {-inf, '-inf', 'min'}
            error('CHEBFUN:SEPARABLEAPPROX:norm:norm', ...
                'SEPARABLEAPPROX does not support this norm.');
            
        case {2, 'op', 'operator'}
            [C, D, R] = cdr( f ); 
            L = C * D * R'; 
            s = svd( L ); 
            normF = s(1);      

        case {'nuc', 'nuclear'}
            [C, D, R] = cdr( f ); 
            L = C * D * R';             
            normF = sum(svd( L )); 
                        
        otherwise
            if ( isnumeric(p) && isreal(p) )
                if ( abs(round(p) - p) < eps )
                    p = round(p); 
                    fp = f.^p;
                    if ( ~mod(p,2) )
                        normF = ( sum2( fp ) ).^( 1/p );
                    else
                        error('CHEBFUN:SEPARABLEAPPROX:norm:norm', ...
                            'p-norm must have p even for now.');
                    end
                else
                    error('CHEBFUN:SEPARABLEAPPROX:norm:norm', ...
                        'SEPARABLEAPPROX does not support this norm.');
                end
            else
                error('CHEBFUN:SEPARABLEAPPROX:norm:unknown', 'Unknown norm.');
            end
            
    end
end

end

%%% 
% Below is the power_method proposed by Mario Bebendorf. This is occasionally
% faster than the current. It is kept here because it may be useful for later.
%%% 
% function [normF, normloc] = power_method( f )
%     % Maximise a CHEBFUN2 using a fast low rank power method.
%     [C, D, R] = cdr( f );
%     n = max( length( C ), length( R ) );
%     x = ones(n,1) ./ n;
%     y = ones(n,1) ./ n;
%     x = x ./ norm( x );
%     y = y ./ norm( y );
%
%     C = C * D; R = R.';
%     rk = length(f);
%     onesrk=ones(1,rk);
%
%     xpts = chebpts(n);
%     C = C(xpts,:); R = R(xpts,:);
%
%     % Perform a fast power iteration.
%     mylam=abs(U(1,1)); xi=1:n; yi=1:n;
%     for jj = 1:1e4
%         C2 = C(xi,:).*x(xi,onesrk);
%         R2 = R(yi,:).*y(yi,onesrk);
%
%         % two norm of low rank vector is frobenius norm of matrix.
%         CC = C2.'*C2; RR = R2.'*R2;
%         myAs = sqrt(sum(sum(CC.*RR)));
%
%         cc = x(xi).'*C2;
%         rr = y(yi).'*R2;
%         mysAs = cc.'*rr; mysAs=mysAs(1,1);
%
%         % two norm of rank-1 vector is Frob norm of rank-1 matrix.
%         xx = x(xi).'*x(xi); yy=y(yi).'*y(yi);
%         myss = xx*yy;
%
%         tmp = mysAs/myss;
%         if abs(mylam-tmp)<1e-5, break; end % stopping criterion.
%         mylam = tmp;
%
%         X = C2/myAs; Y = R2/myAs;
%
%         % Prevent vector increasing in rank by doing the SVD.
%         [Qc, Rc] = qr(X,0); [Qr, Rr] = qr(Y,0);
%         [U, ignored, V] = svd(Rc*Rr.');
%         U = Qc*U; V = Qr*V;
%         x(xi) = U(:,1); y(yi) = V(:,1);
%
%         % update xkeep and ykeep
%         xi=find(abs(x)>100*eps); yi=find(abs(y)>100*eps);
%     end
%
%     % % Estimate maximum from final vectors.
%     [ignored, idx]=max(abs(x)); ystar = xpts(idx);
%     [ignored, idy]=max(abs(y)); xstar = xpts(idy);
%
%     % get domain.
%     rect = f.corners;
%
%     % Do constrained trust region algorithm to finish off.
%     sgn = - sign(feval(f,xstar,ystar));
%
%     try
%         warnstate = warning('off','CHEBFUN2:norm');
%         warning('off')
%         options = optimset('Display','off','TolFun', eps, 'TolX', eps);
%         [Y,m,flag,output] = fmincon(@(x) sgn*feval(f,x(1),x(2)),[xstar,ystar],[],[],[],[],rect([1;3]),rect([2;4]),[],options);
%         warning('on')
%         warning(warnstate)
%     catch
%         % nothing is going to work, so just return the initial guesses.
%         m = abs(feval(f,xstar,ystar));
%         Y = [xstar,ystar];
%         warning('CHEBFUN2:MINANDMAX2','Unable to find Matlab''s optimization toolbox so results will be inaccurate.');
%     end
%     m=-m; normF = m;
%     normloc=Y;
% end
