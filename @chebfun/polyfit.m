function f = polyfit( y, n )  
% POLYFIT Fit polynomial to a chebfun.
%
% F = POLYFIT(Y,N) returns a chebfun F corresponding to the polynomial 
% of degree N that fits the chebfun Y in the least-squares sense.
%
% F = POLYFIT(X,Y,N,D) returns a chebfun F on the domain D which 
% corresponds to the polynomial of degree N that fits the data (X,Y) 
% in the least-squares sense.
%
% Note CHEBFUN/POLYFIT does not not support more than one output argument
% in the way that MATLAB/POLYFIT does.
%
% If y is a global polynomial of degree n then this code has an O(n (log n)^2)
% complexity. If y is piecewise polynomial then it has an O(n^2) complexity.  
%
% See also POLYFIT, DOMAIN/POLYFIT, LEG2CHEB, CHEB2LEG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargout > 1 )
    error('CHEBFUN:polyfit:nargout','Chebfun/polyfit only supports one output');
end

for k = 1 : numel( y )
    f(k) = columnfit(y(k), n);
end

function f = columnfit( y, n )
 
if n > length(y) && y.nfuns == 1
    f = y;
else
    % The code below is a fast version of the code: 
    % Enorm = legpoly(0:n,[a,b],'norm');  % Legendre-Vandermonde matrix   
    % f = Enorm*(Enorm'*y);               % least squares chebfun
      
    [a,b] = domain(y);                              % domain
    if ( y.nfuns > 1 ) 
        scl = 2./(2*(0:n)'+1);                      % orthonormal scaling
        ends = y.ends; yfuns = y.funs;              % piecewise pieces
        cleg = zeros(n+1,1);                        
        for jj = 1:y.nfuns
            yfun = yfuns( jj );                     % For each piece calc
            dom = [ends(jj) ends(jj+1)];            % int P_k f(x)dx, over
            sdom = 2*(dom - a)./(b - a) - 1;        % subdomain. 
            [xscaled, w] = legpts(max(length(yfun),n+1), sdom);
            val = feval( yfun, (xscaled + 1) * (b-a)/2 + a );
            cleg( 1 ) = cleg( 1 ) + w * val;
            cleg( 2 ) = cleg( 2 ) + w * ( xscaled .* val );
            Pm2 = 1; Pm1 = xscaled; 
            for kk = 1:n-1                          % Eval legpoly by rec.
               P = (2-1/(kk+1))*Pm1.*xscaled - kk/(kk+1)*Pm2;  
               cleg( kk+2 ) = cleg( kk+2 ) + w * ( P.*val ); 
               Pm2 = Pm1; Pm1 = P;
            end    
        end
        clegnorm = cleg ./ scl; 
        c = leg2cheb( flipud( clegnorm ) );         % Apply Enorm
        f = chebfun( c, 'coeffs', [a,b] );          % Make a chebfun
    else
        % if y only has one fun, then we can do this much faster. 
        c = y.funs.coeffs;
        c = cheb2leg( c ); 
        c = c( end-n : end ); 
        c = leg2cheb( c ); 
        f = chebfun( c, 'coeffs', [a,b] ); 
    end
end
