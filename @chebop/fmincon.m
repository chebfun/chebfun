function [u,w,gval] = fmincon(g,f,u0,w0,gu,gw) 
% FMINCON 
%    FMINCON attempts to solve problems of the form:
%     min G(U,W)  subject to:  F(U,W) = 0
%     U,W                     
%     
%    See also .

% call gradient descent if no partials are supplied
if ( nargin < 6 )
    [u,w,gval,res] = gradientDescent(g,f,u0,w0);

    [u,w,gval,res] = bfgs(g,f,u,w,res);

% use semi-Newton method if partials are supplied
else
    [u,w,gval,res] = seminewton(g,f,u0,w0,gu,gw);

    [u,w,gval,res] = bfgs(g,f,u,w,res);

end
clf, semilogy(res,'.-')

end



function [u,w,gval,res] = bfgs(g,f,u0,w0,res)
% initialize U and W
w = w0;

% create function handle for local operator
dom = f.domain;
fop = f.op; 
flbc = f.lbc;
frbc = f.rbc;
fbc = f.bc;
N = chebop(dom);

% create local chebop
N.op = @(x,u) fop(x,u,w);
if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

% update u
u = N\0;

% initialize cost
gval = sum(g(u,w));

% initialize bfgs variables
H = @(v) v;
K = 0*w*w';

% do one step of steepest descent without line search
dg = gradient(g,f,u,w);
s = -simplify(H(dg) + K*dg); 
w = simplify(w + s);
  
dgold = dg; 
res = [res,norm(dg,inf)/norm(w,inf)];

% create local chebop
N.op = @(x,u) fop(x,u,w);
if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

% update u
u = N\0;

% initialize cost
gval = sum(g(u,w));

% use adjoint method to compute gradient for bfgs
disp('bfgs')
disp('it      g          grad        delta')
itmax = 200;
tol = 1e-12;
for kk = 1:itmax

    % compute gradient 
    dg = gradient(g,f,u,w);

    % rank two update to inverse Hessian (this might not be most efficient)
    y = simplify(dg - dgold);
    dgold = dg;
    rho = 1/(y'*s);
    hy = H(y) + K*y; 
    K = K + -rho*(s*hy'+hy*s') + (rho+rho*rho*(y'*hy))*s*s';

    % update w
    s = -simplify(H(dg) + K*dg);
    w = simplify(w + s);

    % limit memory of K
    %if ( mod(kk,10) == 0 ), K = 0*K; end

    % create local chebop
    N.op = @(x,u) fop(x,u,w);
    if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
    if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
    if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

    % update u
    u = N\0;

subplot(3,1,1), plot([u,w])
subplot(3,1,2), plot(dg)
subplot(3,1,3), plot(s), pause(.1)

    % update gval
    gval = sum(g(u,w));
disp(sprintf(' %d %e %e %e',kk,gval,norm(dg,inf),norm(s,inf)))

    % update res
    res = [res,norm(dg,inf)/norm(w,inf)];

    % exit if converged
    if ( res(end) < tol ), break, end

end

end





function [u,w,gval,res] = seminewton(g,f,u0,w0,gu,gw)
% initialize U and W
u = u0;
w = w0;
gval = sum(g(u,w));

% create function handle for local operator
dom = f.domain;
fop = f.op; 
flbc = f.lbc;
frbc = f.rbc;
fbc = f.bc;
N = chebop(dom);


% loop for newton iterations
disp('semi-Newton')
disp('it      g          grad        delta')
lam = 0*w;
res = [];
itmax = 20;
tol = 1e-2;
for kk = 1:itmax

    % linearize f
    fpart = linearize(f, [u;w]);
    funs = fpart.constraint.functional;
    vals = fpart.constraint.values;
    fu = linop(fpart{1});
    fu.constraint.functional = funs(:,1);
    fu.constraint.values = 0*vals(:,1);
    fw = linop(fpart{2});
    fw.constraint.functional = funs(:,2);
    fw.constraint.values = 0*vals(:,1);

    % compute adjoint of fu
    [fu_star, fu_star_op, fu_star_bcOpL, fu_star_bcOpR, fu_star_bcOpM] = ...
                                                      adjoint(fu, 'bvp');

    % compute adjoint of fw
    [fw_star, fw_star_op, fw_star_bcOpL, fw_star_bcOpR, fw_star_bcOpM] = ...
                                                      adjoint(fw, 'bvp');

    % create local chebop
    N.op = @(x,u,w,lam) [ gu(u,w) + fu_star_op(x,lam);...
                          gw(u,w) + fw_star_op(x,lam);...
                          fop(x,u,w) ];
    
    % left boundary conditions
    lbcs = [];
    if ( ~isempty(fu_star_bcOpL) && ~isempty(lbcs) )
       lbcs = @(u,w,lam) [lbcs(u,w,lam);fu_star_bcOpL(lam)]; 
    elseif ( ~isempty(fu_star_bcOpL) && isempty(lbcs) )
       lbcs = @(u,w,lam) [fu_star_bcOpL(lam)]; 
    end
    if ( ~isempty(fw_star_bcOpL) && ~isempty(lbcs) )
       lbcs = @(u,w,lam) [lbcs(u,w,lam);fw_star_bcOpL(lam)]; 
    elseif ( ~isempty(fw_star_bcOpL) && isempty(lbcs) )
       lbcs = @(u,w,lam) [fw_star_bcOpL(lam)]; 
    end
    if ( ~isempty(flbc) && ~isempty(lbcs) )
       lbcs = @(u,w,lam) [lbcs(u,w,lam);flbc(u,w)]; 
    elseif ( ~isempty(flbc) && isempty(lbcs) )
       lbcs = @(u,w,lam) [flbc(u,w)]; 
    end
    if ( ~isempty(lbcs) )
       N.lbc = @(u,w,lam) lbcs(u,w,lam); 
    end

    % right boundary conditions
    rbcs = [];
    if ( ~isempty(fu_star_bcOpR) && ~isempty(rbcs) )
       rbcs = @(u,w,lam) [rbcs(u,w,lam);fu_star_bcOpR(lam)]; 
    elseif ( ~isempty(fu_star_bcOpR) && isempty(rbcs) )
       rbcs = @(u,w,lam) [fu_star_bcOpR(lam)]; 
    end
    if ( ~isempty(fw_star_bcOpR) && ~isempty(rbcs) )
       rbcs = @(u,w,lam) [rbcs(u,w,lam);fw_star_bcOpR(lam)]; 
    elseif ( ~isempty(fw_star_bcOpR) && isempty(rbcs) )
       rbcs = @(u,w,lam) [fw_star_bcOpR(lam)]; 
    end
    if ( ~isempty(frbc) && ~isempty(rbcs) )
       rbcs = @(u,w,lam) [rbcs(u,w,lam);frbc(u,w)]; 
    elseif ( ~isempty(frbc) && isempty(rbcs) )
       rbcs = @(u,w,lam) [frbc(u,w)]; 
    end
    if ( ~isempty(rbcs) )
       N.rbc = @(u,w,lam) rbcs(u,w,lam); 
    end

    % mixed boundary conditions
    bcs = [];
    if ( ~isempty(fu_star_bcOpM) && ~isempty(bcs) )
       bcs = @(x,u,w,lam) [bcs(x,u,w,lam);fu_star_bcOpM(x,lam)]; 
    elseif ( ~isempty(fu_star_bcOpM) && isempty(bcs) )
       bcs = @(x,u,w,lam) [fu_star_bcOpM(x,lam)]; 
    end
    if ( ~isempty(fw_star_bcOpM) && ~isempty(bcs) )
       bcs = @(x,u,w,lam) [bcs(x,u,w,lam);fw_star_bcOpM(x,lam)]; 
    elseif ( ~isempty(fw_star_bcOpM) && isempty(bcs) )
       bcs = @(x,u,w,lam) [fw_star_bcOpM(x,lam)]; 
    end
    if ( ~isempty(fbc) && ~isempty(bcs) )
       bcs = @(x,u,w,lam) [bcs(x,u,w,lam);fbc(x,u,w)]; 
    elseif ( ~isempty(fbc) && isempty(bcs) )
       bcs = @(x,u,w,lam) [fbc(x,u,w)]; 
    end
    if ( ~isempty(bcs) )
       N.rbc = @(x,u,w,lam) bcs(x,u,w,lam); 
    end

    % construct rhs
    rhs = N*[u;w;lam];

    % update res
    res = [res,norm(rhs,inf)];

    % linearize N
    y = [u;w;lam];
    Ny = linearize(N, y);

    % convert chebfun blocks into multiplication operators
    dims = size(Ny);
    for ii = 1:dims(1)
      for jj = 1:dims(2) 
        b = Ny.blocks{ii,jj};
        if ( strcmp(class(b),'chebfun') )
          Ny.blocks{ii,jj} = operatorBlock.mult(b,dom);
        end
      end
    end

    % convert scalars into zero functionals
    funs = Ny.constraint.functional;
    dims = size(funs);
    for ii = 1:dims(1)
      for jj = 1:dims(2) 
        b = funs(ii,jj);
        if ( strcmp(class(b{1}),'double') && b{1} == 0 )
          funs(ii,jj) = functionalBlock.feval(dom(1),dom)*operatorBlock.mult(0,dom);
        elseif ( strcmp(class(b{1}),'double') && b{1} ~= 0 )
          error('Not allowed!')
        end
      end
    end
    Ny.constraint.functional = funs;

    % single newton step
    prefs = cheboppref();
    prefs.discretization = @chebcolloc2;
    dy = linsolve(Ny, rhs, prefs);

subplot(3,1,1), plot(y)
subplot(3,1,2), plot(rhs)
subplot(3,1,3), plot(dy), pause(.1)

    % update u, w, and lam
    u = u - dy{1};
    w = w - dy{2};
    lam = lam - dy{3};

    % update gval
    gval = sum(g(u,w));
disp(sprintf(' %d %e %e %e',kk,gval,norm(rhs,inf),norm(dy,inf)))

    % check for convergence
    if ( norm(rhs,inf) < tol*abs(gval) && kk > 1), break, end

end

end





function [u,w,gval,res] = gradientDescent(g,f,u0,w0)
% initialize U and W
u = u0;
w = w0;

% create function handle for local operator
dom = f.domain;
fop = f.op; 
flbc = f.lbc;
frbc = f.rbc;
fbc = f.bc;
N = chebop(dom);

    % create local chebop
    N.op = @(x,u) fop(x,u,w);
    if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
    if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
    if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

    % update u
    u = N\0;

gval = sum(g(u,w));
scl_init = .1;
scl_max = .5;

% use adjoint method to compute gradient
disp('gradient descent')
disp('it      g          grad        delta')
tol = 1e-2;
itmax = 40;
res = [];
for kk = 1:itmax

    % compute gradient
    dg = gradient(g,f,u,w);

    % update w and u
    if ( kk == 1 )
      w = w - scl_init*dg;
    elseif ( res(end) > scl_init )
      w = w - scl_init*dg;
    else
      w = w - scl_max*dg;
    end

    % create local chebop
    N.op = @(x,u) fop(x,u,w);
    if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
    if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
    if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

    % update u
    u = N\0;
   
    % update gval
    gval = sum(g(u,w));

    % update res
    res = [res,norm(dg,inf)/norm(w,inf)];
disp(sprintf(' %d %e %e %e',kk,gval,norm(dg,inf),res(end)))

subplot(2,1,1),plot([u,w,dg])
subplot(2,1,2),semilogy((1:length(res))',res,'.-'), pause(.1)

    % eit if converged
    if ( res(end) < tol ), break, end

end

end





function dg = gradient(g,f,u,w)

% create function handle for local operator
dom = f.domain;
fop = f.op; 
flbc = f.lbc;
frbc = f.rbc;
fbc = f.bc;
N = chebop(dom);

% use adjoint method to compute gradient 
one = chebfun('1',dom);

    % create local chebop
    N.op = @(x,u) fop(x,u,w);
    if ( ~isempty(flbc) ), N.lbc = @(u) flbc(u,w); end
    if ( ~isempty(frbc) ), N.rbc = @(u) frbc(u,w); end
    if ( ~isempty(fbc) ), N.bc = @(x,u) fbc(x,u,w); end

    % update u
    u = N\0;

    % linearize g
    y = [u;w];
    gpart = linearize(chebop(@(x,u,w) g(u,w), dom), y);
    gu = linop(gpart{1});
    gw = linop(gpart{2});

    % linearize f
    fpart = linearize(f, y);
    funs = fpart.constraint.functional;
    vals = fpart.constraint.values;
    fu = linop(fpart{1});
    fu.constraint.functional = funs(:,1);
    fu.constraint.values = 0*vals(:,1);
    fw = linop(fpart{2});
    fw.constraint.functional = funs(:,2);
    fw.constraint.values = 0*vals(:,1);

    % compute adjoint of fu
    fu_star = adjoint(fu,'bvp');

    % solve for lambda
    prefs = cheboppref();
    prefs.discretization = @chebcolloc2;
    gu = gu(one);
    lam = linsolve(fu_star, gu, prefs);

    % compute gradient 
    fw = fw*1;
    dg = -fw{1}.*lam{1} + gw(one);
    dg = simplify(dg{1});

end
