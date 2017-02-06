function [r, pol, res, zer, z, f, w, errvec, p, q] = aaamn_lawson(Z,F,m,n,tol,Lawsoniter,doplot)
% [r, pol, res, zer, z, f, w, errvec, p, q] = aaamn_lawson(Z,F,m,n,tol,Lawsoniter,doplot)
% near-best rational approximation of data F on set Z of type (MMAX,NMAX). 
% It first finds a AAA rational approximant to F \approx r(Z), then
% attempts to refine the approximant by a Lawson process, i.e. an iterative
% reweighting. 
% 
% This code is originally designed for computing good reference points for the
% rational remez code to follow, but can be used independently for
% constructing a rational approximation r that can be much closer than AAA
% to the best rational approximant. 
%
% Input:  Z = vector of sample points
%         F = vector of data values, or a function handle
%         m, n: max type is (m,n), set to 10 if omitted
%         tol = relative tolerance tol, default: 1e-13 
%         Lawsoniter: max. iteration number of Lawson updates (default 10)
%         doplot: 1 to plot error curve history (default 0)
%         
%
% Output: r = AAA-Lawson approximant to F (function handle)
%         pol,res,zer = vectors of poles, residues, zeros
%         errvec = vector of errors at each step
%         z,f,w = vectors of support pts, function values, weights
%         p,q = chebfuns s.t. r= p/q (note this can be numerically unstable)

M = length(Z);                            % number of sample points
if ( (nargin < 3) || isempty(m) ), m = 10;% default max type (10,10) 
    disp('No input type given; will compute type (10,10) approximant')
end                 
if ( (nargin < 4) || isempty(n) ), n = m; end
if ( (nargin < 5) || isempty(tol) ), tol = 1e-13; end % default relative tol 
mmax = m+1; nmax = n+1;                   % for coding convenience
if ( (nargin < 6) || isempty(Lawsoniter) )% number of Lawson updates
    Lawsoniter = max([5 min([20,mmax,nmax])]); 
end 
if ( (nargin < 7) || isempty(doplot)), doplot = 0;  end  % plot Lawson updates
if ~isfloat(F), F = F(Z); end             % convert function handle to vector
Z = Z(:); F = F(:);                       % work with column vectors
SF = spdiags(F,0,M,M);                    % left scaling matrix
J = 1:M;                                  % indices that are not support pts
z = []; f = []; C = [];                   % initializations
errvec = []; R = mean(F); 
for mn = 1:max(mmax,nmax)
  [~,j] = max(abs(F-R));                  % select next support point
  z = [z; Z(j)];                          % update set of support points
  f = [f; F(j)];                          % update set of data values
  J(J==j) = [];                           % update index vector
  C = [C 1./(Z-Z(j))];                    % next column of Cauchy matrix
  Sf = diag(f);                           % right scaling matrix
  A = SF*C - C*Sf;                        % Loewner matrix
    
           if mn > min(nmax,mmax) % nondiagonal case, find projection subspace 
                if mmax < nmax
                q = f(:);
                else
                q = ones(length(z),1);
                end
                Q = orthspace(z,mn-min(mmax,nmax),q);     % projection subspace                              
                [~,~,V] = svd(A(J,:)*Q,0);              % SVD on projected subspace
                w = Q*V(:,end);
           else             
                [~,~,V] = svd(A(J,:),0);               % SVD, no projection needed
                w = V(:,mn);                           % weight vector             
           end     
  wf = w.*f;
  N = C*(w.*f); D = C*w;                  % numerator and denominator
  R = F; R(J) = N(J)./D(J);               % rational approximation
  err = norm(F-R,inf);
  errvec = [errvec; err];                 % max error at sample points
  if err < tol*norm(F,inf), break, end    % stop if converged
end
r = @(zz) feval(@rr,zz,z,w,f);            % AAA approximant as function handle

Rori = R;
% now start Lawson, in this mode we leave interpolation and work with
% 'alpha-beta' mode. 
          wei = ones(length(J),1);
          nrmbest = inf;
          if mn > min(nmax,mmax)  % Deal with projection for m neq n                                             
              if mn>nmax
                A =[SF*C*Q -C];        
              else % need to redefine Q as not the same as AAA above        
                q = ones(length(z),1);
                Q = orthspace(z,mn-min(mmax,nmax),q);     % projection subspace                              
                A =[SF*C -C*Q];                   
              end
          else
            A =[SF*C -C]; % diagonal case
          end
          
	      for it = 1:Lawsoniter
              weiold = wei; 
              wei = wei .* sqrt(abs(F(J)-R(J)));       % update Lawson weights
              wei = wei/sum(wei);                      % normalize 
              if norm(weiold-wei)/norm(wei)< sqrt(tol) % declare Lawson converged
                  break
              end
              D = spdiags(sqrt(wei),0,length(wei),length(wei)); % diagonal weight matrix

              [~,~,V] = svd(D*A(J,:),0);     % weighted least-squares via SVD
          if mn > min(nmax,mmax)
              if mn > nmax
                w = Q*V(1:nmax,end); wf = V(nmax+1:end,end);            
              else
                w = V(1:nmax,end); wf = Q*V(nmax+1:end,end); 
              end
          else
            w = V(1:mn,end); wf = V(mn+1:2*mn,end);            
          end                      
            f = wf./w;                          % for compatibility with interpolatory-aaa                        
            N = C*wf; D = C*w;                  % numerator and denominator               
            R = F; R(J) = N(J)./D(J);           % rational approximation
            err = norm(F-R,inf);            
            errvec = [errvec; err];                 % max error at sample points                            
            if err < nrmbest        % adopt best so far
            nrmbest = norm(F-R,'inf'); 
            r = @(zz) feval(@rrab,zz,z,w,wf,f); % AAA approximant as function handle                      
            end
                        
            if doplot
                subplot(2,1,1)
                plot(Z,F-Rori,'r.','markersize',12,'linewidth',.5)
                title('AAA error')
                grid on, hold on
                if exist('hh','var'), set(hh,'color',(0.8)*[1 1 1]); end
                subplot(2,1,2)
                title('AAA-Lawson error')
                hh = plot(Z,F-R,'k.','markersize',12,'linewidth',.5);                        
                grid on, hold on
                plot(z,0*z,'m.','markersize',14)
                ylim(err*[-1 1]); drawnow, shg            
            end
            
          end          

    % compute poles and roots
    B = eye(mn+1); B(1,1) = 0;                 
    E = [0 w.'; ones(mn,1) diag(z)];      
    zer = eig(E,B); zer = zer(~isinf(zer));   % zeros  
    E = [0 w.'; ones(mn,1) diag(z)];      
    pol = eig(E,B); pol = pol(~isinf(pol));   % poles
    dz = 1e-5*exp(2i*pi*(1:4)/4);
    res = r(bsxfun(@plus,pol,dz))*dz.'/4;     % residues
    
    if ( nargout > 8 )  % form p and q via N/D, NOTE: not always stable
    D = @(x) 0;    N = @(x) 0;     % form r=N/D in barycentric form
    for ii = 1:length(z)
       D = @(x) D(x) + w(ii)./(x-z(ii));
       N = @(x) N(x) + wf(ii)./(x-z(ii));               
    end
    dom = [min(real(Z)) max(real(Z))];  % set domain
    x = chebpts(mmax+nmax+1,dom);
    node = @(x) prod(x-z);              % needed to check sign
    nodex = zeros(length(x),1);         % setup node values    
    for ii = 1:length(x),    nodex(ii) = node(x(ii)); end
    qvals = nodex.*feval(D,x);          % values of p,q
    pvals = nodex.*feval(N,x);
    p = chebfun(pvals,dom); q = chebfun(qvals,dom); % form chebfuns    
    end
end    



function r = rr(zz,z,w,f)                 % evaluate r at zz
zv = zz(:);                               % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.');            % Cauchy matrix 
r = (CC*(w.*f))./(CC*w);                  % AAA approx as vector
ii = find(isnan(r));                      % Find values NaN = 0/0 if any
for j = 1:length(ii)
  r(ii(j)) = f(find(zv(ii(j))==z));       % Force interpolation there
end
r = reshape(r,size(zz));                  % AAA approx
end

function r = rrab(zz,z,w,wf,f)                 % evaluate r at zz
zv = zz(:);                               % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.');            % Cauchy matrix 
r = (CC*(wf))./(CC*w);                  % AAA approx as vector

ii = find(isnan(r));                      % Find values NaN = 0/0 if any
for j = 1:length(ii)
  r(ii(j)) = f(find(zv(ii(j))==z));       % Force interpolation there
end

r = reshape(r,size(zz));                  % AAA approx
end

function Q = orthspace(z,dim,q)       % orthonormal projection space for (m,n)
if ( dim == 0 ), Q = eye(length(z)); end 
if ( nargin < 3 ), q = ones(length(z),1); end
                Q = q/norm(q);
                for ii = 2:dim               % orthogonal complement via Lanczos-type process
                Qtmp = diag(z)*Q(:,end);
                Qtmp = Qtmp - Q*(Q'*Qtmp);   % orthogonalize
                Qtmp = Qtmp/norm(Qtmp);      % normalize
                Q = [Q Qtmp]; 
                end
            [Q,~] = qr(Q); Q = conj(Q(:,dim+1:end)); % desired null space
end