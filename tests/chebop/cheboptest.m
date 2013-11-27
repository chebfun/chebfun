function [u, info] = cheboptest(problemno)
cheboppref('display','iter')
cheboppref('plotting','on')
cheboppref('damped','on')
switch problemno
    case 1
        N = chebop(@(x,u) diff(u,2) + sin(u));
        N.lbc = @(u) u-2; N.rbc = @(u) u -2;
        N.init = chebfun(@(x) 0*x+2);
        rhs = 0;
    case 2
        cheboppref('damped','on')
        N = chebop(@(x,u) 0.05*diff(u,2) + u.^2 - 1);
        N.lbc = @(u) u; N.rbc = @(u) u - 1;
        N.init = chebfun(@(x) .5*(x+1));
        rhs = 0;
    case 3
        N = chebop(@(x,u) 0.01*diff(u,2) + 2*(1-x.^2).*u + u.^2);
        N.lbc = @(u) u;
        N.rbc = @(u) u;
        x = chebfun('x');
        rhs = 1;
        N.init = 0*2*(x.^2-1).*(1-2./(1+20*x.^2));
    case 4
        dom = [0 1];
        N = chebop(@(x,u) diff(u,2)-8.*sinh(8.*u), dom);
        N.lbc = @(u) u; N.rbc = @(u) u - 1;
        N.init = chebfun(@(x) x, dom);
        rhs = 0;
    case 5
        % Problem with a breakpoint
        N = chebop(@(x,u) diff(u,2) + sign(x).*sin(u));
        N.lbc = @(u) u-2; N.rbc = @(u) u -2;
        N.init = chebfun(@(x) 0*x+2);
        rhs = 0;
        
end
[u, info] = solvebvp(N, rhs);