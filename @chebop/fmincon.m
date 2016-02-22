function [u,w,gval] = fmincon(G,F,u0,w0) 
% FMINCON 
%    FMINCON attempts to solve problems of the form:
%     min G(U,W)  subject to:  F(U,W) = 0
%     U,W                     
%     
%    See also .

% initialize U and W
u = u0;
w = w0;
gval = sum(G*[u;w]);
scl = .01;

% create function handle for local operator
Fop = F.op; 
Flbc = F.lbc;
Frbc = F.rbc;
Fbc = F.bc;
N = chebop(F.domain);

% use adjoint method to compute gradient
one = chebfun('1',F.domain);
itmax = 200;
for ii = 1:itmax

    % create local chebop
    N.op = @(x,u) Fop(x,u,w);
    if ( ~isempty(Flbc) ), N.lbc = @(u) Flbc(u,w); end
    if ( ~isempty(Frbc) ), N.rbc = @(u) Frbc(u,w); end
    if ( ~isempty(Fbc) ), N.bc = @(x,u) Fbc(x,u,w); end

    % linearize G
    y = [u;w];
    Gpart = linearize(G, y);
    Gu = linop(Gpart{1});
    Gw = linop(Gpart{2});

    % linearize F
    Fpart = linearize(F, y);
    funs = Fpart.constraint.functional;
    vals = Fpart.constraint.values;
    Fu = linop(Fpart{1});
    Fu.constraint.functional = funs(:,1);
    Fu.constraint.values = 0*vals(:,1);
    Fw = linop(Fpart{2});
    Fw.constraint.functional = funs(:,2);
    Fw.constraint.values = 0*vals(:,1);

    % compute adjoint of Fu
    Fu_star = adjoint(Fu,'bvp');

    % update u
    u = N\0;

    % solve for lambda
    prefs = cheboppref();
    prefs.discretization = @chebcolloc2;
    gu = Gu(one);
    lam = linsolve(Fu_star, gu, prefs);

    % compute gradient 
    fw = Fw*1;
    dG = -fw{1}.*lam{1} + Gw(one);
    dG = simplify(dG{1});

    % update w and u
    w = w - scl*dG;
    gval = sum(G*[u;w]);

[gval,norm(dG,inf)]
plot([u,w,dG]), pause(.1)

end


end
