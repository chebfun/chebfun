function [r, pol, res, zer, z, Z, f, w, wf, errvec, p, q] = aaamn_lawson(F, varargin)
%AAAMN_Lawson   near-best rational approximation of F. 
% 
% R = aaamn_lawson(F) computes a rational aproximant of (default) type
% (10,10) on the default interval [-1,1]
%
% [R, POL, RES, ZER] = aaamn_lawson(...) outputs the poles, residues and zeros
% of R. The poles, zeros will approximate those of F (not well if R-F is not small)
% 
% [R, POL, RES, ZER, z, Z, F, W,WF, ERRVEC, P, Q] = aaamn_lawson(...) 
% outputs additionally the sample points Z, support points z, values
% F=f(Z), weights w and wf (see below) and AAA errvec, and p,q = chebfuns 
% s.t. r= p/q (note this can be numerically unstable)
%
% [...] = aaamn_lawson(F,m,n) specifies the type to (m,n).
%
% [...] = aaamn_lawson(F,Z,m,n) also specifies the sample points Z
% (recommended). 
%
% [...] = aaamn_lawson(F,Z,m,n,'plot','on') will plot the error functions
% as the Lawson iterations proceed. 
%
% [...] = aaamn_lawson(F,m,n,'dom',[-1,2]) specifies the domain (this has no effect 
% if Z is specified)
%
% [...] = aaamn_lawson(F,m,n,'tol',1e-5) specifies the Lawson iterate
% relative stopping criterion (here to 1e-5)
%
% [...] = aaamn_lawson(F,m,n,'iter',10) limits the Lawson maximum
% iterations to 10. 
% 
% 
% The algorithm first finds a AAA rational approximant to F \approx r(Z), then
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
%     R = AAAMN_Lawson(F, Z, m, n, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance for Lawson iteration (default 1e-5)
%   - 'iter', IT: maximal number of Lawson iterations (default MMAX = max([5 min([20, m, n])])).
%         
%
% Output: r = AAA-Lawson approximant to F (function handle)
%         pol,res,zer = vectors of poles, residues, zeros
%         errvec = vector of errors at each step
%         z,f,w,wf = vectors of support pts, function values, weights
%         s.t. r = N/D, N(x) = sum_i wf(i)/(x-z(i)) and 
%         D(x) = sum_i w(i)/(x-z(i)).
%         p,q = chebfuns s.t. r= p/q (note this can be numerically unstable)
%
% Examples:
%    r = aaamn_lawson(@abs,10,10)
%    r = aaamn_lawson(@abs,10,10,'plot','on')
%    [r, pol, res, zer, z, f, w, wf, errvec, p, q] = aaamn_lawson(@abs,10,10,'plot','on','dom',[-1 2])
%    r = aaamn_lawson(@exp,4,2,'plot','on','dom',[-1 2])
%    r = aaamn_lawson(@(x)log(1.1-x),5,5,'plot','on')
%
%    f = chebfun(@(x)-1./(log(abs(x)).^2),[-.1,.1],'splitting','on'); 
%    [r,pol,res] = aaamn_lawson(f,linspace(-.1,.1,1e4),18,18,'plot','on')
%
%
%
%   Reference:
%   [1] Yuji Nakatsukasa, Olivier Sete, Lloyd N. Trefethen, "The AAA algorithm
%   for rational approximation", arXiv:1612.00337.
%
%   [2] Bernhard Beckermann, Silviu Filip, Yuji Nakatsukasa, rational
%   remez paper, in preparation (2017).  
%
% See also AAA, CF, CHEBPADE, PADEAPPROX, RATINTERP, REMEZ

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% parse inputs
[F, Z, m, n, Lawsoniter, tolLawson, doplot, tol ] = ...
    parseInputs(F, varargin{:});

M = length(Z);                            % number of sample points
mmax = m+1; nmax = n+1;                   % for coding convenience
if ( (nargin < 6) || isempty(Lawsoniter) )% number of Lawson updates
    Lawsoniter = max([5 min([20,mmax,nmax])]); 
end 
if ~isfloat(F), F = feval(F,Z); end        % convert function handle to vector
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
    
        if ( mn > min(nmax,mmax) ) % nondiagonal case, find projection subspace 
             if mmax < nmax
             q = f(:);
             else
             q = ones(length(z),1);
             end
             Q = orthspace(z,mn-min(mmax,nmax),q);   % projection subspace 
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
  if ( err < tol*norm(F,inf) ), break, end    % stop if converged
end
    r = @(zz) feval(@rr,zz,z,w,f);            % AAA approximant as function handle
    Rori = R;
    
% now start Lawson, in this mode we leave interpolation and work with
% 'alpha-beta' mode. 
          wei = ones(length(J),1);
          nrmbest = inf;
          if ( mn > min(nmax,mmax) )  % Deal with projection for m neq n                                             
              if mn>nmax
                A =[SF*C*Q -C];        
              else % need to redefine Q as not the same as AAA above        
                q = ones(length(z),1);
                Q = orthspace(z,mn-min(mmax,nmax),q);     % projection subspace                              
                A =[SF*C -C*Q];                   
              end
          else
            A =[SF*C -C];   % diagonal case
          end
          
          rate = 1;         % default Lawson rate, will shrink if not converging
          nrmincreased = 0; % initialization    
	      for it = 1:Lawsoniter
              weiold = wei; 
              wei = wei .* power(abs(F(J)-R(J)),rate); % update Lawson weights
              wei = wei/sum(wei);                      % normalize 
              if norm(weiold-wei)/norm(wei)< tolLawson % declare Lawson converged
                  break
              end
              D = spdiags(sqrt(wei),0,length(wei),length(wei)); % diagonal weight matrix

              [~,~,V] = svd(D*A(J,:),0);     % weighted least-squares via SVD
              
          if ( mn > min(nmax,mmax) )    % deal with nondiagonal case
              if ( mn > nmax )
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
            errvec = [errvec; err];             % max error at sample points                            
            if ( err < nrmbest )    % adopt best so far
            nrmbest = norm(F-R,'inf'); 
            weibest = wei;                      % store best weight
            r = @(zz) feval(@rrab,zz,z,w,wf,f); % AAA approximant as function handle
            else
                nrmincreased = nrmincreased + 1;
            end
            if ( nrmincreased >= 3 )     % perhaps not converging,
            rate = max( rate/2,0.01 );   % make Lawson update conservative
            if doplot
                warning(['Lawson rate made conservative to ',num2str(rate)])
            end
              nrmincreased = 0; 
            end

            if doplot  % plot error functions (hopefully near-equioscillating)
                subplot(2,1,1)
                plot(Z,F-Rori,'r.','markersize',8)
                title('AAA error')
                grid on, hold on
                if exist('hh','var')
                    set(hh,'color',(0.8)*[1 1 1]); 
                end
                subplot(2,1,2)
                title('AAA-Lawson error')
                hh = plot(Z,F-R,'k.','markersize',8);                        
                grid on, hold on
                h2 = plot(z,0*z,'m.','markersize',12);
                ylim(err*[-1 1]); drawnow, shg            
                if ( it == Lawsoniter ) % plot best function 
                plot(Z,F-r(Z),'b.','markersize',10);
                end
            end            
	      end          

    % compute poles and roots
    B = eye(mn+1); B(1,1) = 0;                 
    E = [0 wf.'; ones(mn,1) diag(z)];      
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
    for ii = 1:length(x), nodex(ii) = node(x(ii)); end
    qvals = nodex.*feval(D,x);          % values of p,q
    pvals = nodex.*feval(N,x);
    p = chebfun(pvals,dom); q = chebfun(qvals,dom); % form chebfuns    
    end
end    


%% parse Inputs:

function [F, Z, m, n, Lawsoniter, tolLawson, doplot, tol ] = ...
    parseInputs(F, varargin)
% Input parsing for AAAmn_lawson.

% Check if F is empty:
if ( isempty(F) )
    error('CHEBFUN:aaamn_lawson:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('CHEBFUN:aaamn_lawson:nColF', 'Input chebfun must have one column.')
    end
end

% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    if length(varargin{1})>2   % sample points Z given. 
    Z = varargin{1};
    varargin(1) = [];    
    else                       % sample points not given.
        
    end
end

% m,n
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    if length(varargin{1})>1 % input (f,Z,[m n])
    mn = varargin{1};
    m = mn(1); n = mn(2);
    varargin(1) = [];    
    elseif isfloat(varargin{2}) % input (f,Z,m,n)
    m = varargin{1};        
    n = varargin{2};    
    varargin([1, 2]) = [];
    else
    m = varargin{1};        
    n = m; % input (f,Z,m), default to diagonal type (m,m)
    varargin(1) = [];
    end
end

if ( ~exist('m', 'var') ) 
     warning('CHEBFUN:aaamn_lawson: type (m,n) not specified, default to (10,10)')
     m = 10; n = 10; 
end

% Set defaults for other parameters:
tolLawson = 1e-5;                       % Relative tolerance for Lawson update.
tol = 1e-15;                            % AAA tolerance
Lawsoniter = max([5 min([20, m, n])]);  % Maximum number of terms.
doplot = 0;                             % Don't plot intermediate functions unless specified

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) |  strncmpi(varargin{1}, 'tolLawson', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tolLawson = varargin{2};   % Lawson tolerance
        end
        varargin([1, 2]) = [];

    elseif ( strncmpi(varargin{1}, 'iter', 4) | strncmpi(varargin{1}, 'maxit', 5))
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            Lawsoniter = varargin{2};  % maximum Lawson iterations
        else
        warning(['CHEBFUN:aaamn_lawson:iter unspecified, use default itermax ', num2str(Lawsoniter)])
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('CHEBFUN:aaamn_lawson:dom', ...
                    ['Given domain does not match the domain of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif strncmpi(varargin{1}, 'plot', 4)  % plot error functions
        if isfloat(varargin{2})
            doplot = varargin{2};
        elseif ( strncmpi(varargin{2}, 'true', 4) | strncmpi(varargin{2}, 'on', 2) )
            doplot = 1;
        end
        varargin([1, 2]) = [];                
    else
        error('CHEBFUN:aaamn_lawson:UnknownArg', 'Argument unknown.')
    end
end

% Deal with Z and F:
if ( ~exist('Z', 'var') && isfloat(F) )
    % F is given as data values, pick same number of sample points:
    Z = linspace(dom(1), dom(2), length(F)).';
end

if ( exist('Z', 'var') )
    % Work with column vector:
    Z = Z(:);
    M = length(Z);
    
    % Function values:
    if ( isa(F, 'function_handle') || isa(F, 'chebfun') )
        % Sample F on Z:
        F = F(Z);
    elseif ( isnumeric(F) )
        % Work with column vector and check that it has correct length.
        F = F(:);
        if ( length(F) ~= M )
            error('CHEBFUN:aaamn_lawson:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    else
        error('CHEBFUN:aaamn_lawson:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    % in AAA this is done adaptively. This can be done with Lawson, but
    % probably safe to take as many points as reasonably possible here. 
    Z = linspace(dom(1), dom(end), 4000).';    
end

end % End of PARSEINPUT().


%% generate function handle, interpolatory mode
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

%% generate function handle, non-interpolatory mode
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

%% null space for m~=n
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

