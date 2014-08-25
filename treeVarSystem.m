dom = [0 20];
t = chebfun(@(t) t, dom);
u = treeVar([1 0 0]);
v = treeVar([0 1 0]);
w = treeVar([0 0 1]);
initCon = [-14 -15 20];
resFun = @(t, u, v, w) [diff(u) + 10*u - 10*v, diff(v) + v - 28*u + u.*w, ...
    diff(w) + (8/3).*w - u.*v]; 
resSys = resFun(t, u, v, w);

systemInfix = cell(length(resSys),1);
coeffs = systemInfix;
varArrays = systemInfix;

numArgs = length(resSys(1).tree.diffOrder);

treeVar.plotTree(resSys(1).tree)

% First look at all diffOrders to ensure we start with the correct indices
indexStart = zeros(1, numArgs);
indexStartDer = indexStart;
totalDiffOrders = indexStart;
for wCounter = 1:length(resSys)
    newIndex =  [1 (cumsum(resSys(wCounter).tree.diffOrder(1:end-1)) + (1:(numArgs-1)))];
    newIndexDer =  ...
        [1 (cumsum(resSys(wCounter).tree.diffOrder(1:end-1) + 1) + (1:(numArgs-1)))];
    indexStart = max(indexStart, newIndex);
    indexStartDer = max(indexStartDer, newIndexDer);
    totalDiffOrders = max(totalDiffOrders, resSys(wCounter).tree.diffOrder);
end
coeffArg = zeros(1, indexStartDer(end) + totalDiffOrders(end));
%%
t = chebfun(@(t) t);
for wCounter = 1:length(resSys)
    
    res = resSys(wCounter);
    diffOrders = res.tree.diffOrder;
    
    expTree = treeVar.expandTree(res.tree, diffOrders);
    
    [newTree, derTree] = treeVar.splitTree(expTree, diffOrders);
    
    [infixDer, dummy, varArrayDer] = treeVar.tree2infix(derTree, wCounter, indexStartDer);
    
    % Find what argument corresponds to the highest derivative one:
    maxDerLoc = find(expTree.diffOrder == max(diffOrders));
    coeffFun = treeVar.toAnon(infixDer, varArrayDer);
    
    % Reset coeffArg:
    coeffArg = 0*coeffArg;
    % Replace one of the 0 in coeffFun with 1 so that we can evaluate COEFFFUN:
    if ( maxDerLoc == numArgs )
        coeffArg(end) = 1;
    else
        % The variable with index maxDerLoc+1 is the next variable we need to
        % start numbering at. So subtract 1 for the index of the highest
        % derivate we're currently interested in.
        coeffArg(indexStartDer(maxDerLoc+1) - 1) = 1;
    end
    coeffs{wCounter} = coeffFun(t, coeffArg);
    % Need to negate the syntax tree as we're moving it to the right-hand side. But
    % if it already starts with a unary minus, we can simply remove it rather than
    % doing a double negation:
    newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
    [infix, varCounter, varArray] = treeVar.tree2infix(newTree, wCounter, indexStart);
    systemInfix{wCounter} = infix;
    varArrays{wCounter} = varArray;
end
%%
coeffs
% coeffs{1} = 1; coeffs{2} = 1; coeffs{3} = 1;
funOut = treeVar.toRHS(systemInfix, varArrays, coeffs, indexStart, totalDiffOrders);
opts = odeset('abstol',1e-13,'reltol',1e-13);
sol = ode113(funOut,[0 10], initCon, opts)
x = sol.x;
y = sol.y;
plot(x, y)

%%
figure
LW = 'linewidth'; FS = 'fontsize';
plot3(y(1,:),y(2,:),y(3,:), LW, 1.6), view(20,20)
axis([-20 20 -40 40 5 45]), grid on
xlabel 'x(t)', ylabel 'y(t)', zlabel 'z(t)'
title('A 3D Trajectory of the Lorenz Attractor - Standard ODE113', FS, 14)

%% With chebfun
figure
u = chebfun.ode113(funOut,[0,10], initCon, opts)
LW = 'linewidth'; FS = 'fontsize';
plot3(u(:,1),u(:,2),u(:,3), LW, 1.6), view(20,20)
axis([-20 20 -40 40 5 45]), grid on
xlabel 'x(t)', ylabel 'y(t)', zlabel 'z(t)'
title('A 3D Trajectory of the Lorenz Attractor - Chebfun solution', FS, 14)

%% With CHEBOP -- Lotka-Volterra
dom = [0 15];
N = chebop(@(t,u,v) [diff(u) - u + u.*v; diff(v) + v - u.*v], dom);
N.lbc = @(u,v) [v-1.2; u-1.2];
uv = N\[0;0]
plot(uv)

%% With CHEBOP -- Higher order Lotka-Volterra, useful for checking BC parsing
dom = [0 3];
N = chebop(@(t,u,v) [diff(u, 2) - u + u.*v; diff(v) + v - u.*v], dom);
N.lbc = @(u,v) [u-1; v-1.5; diff(u)-0.5];
uv = N\[0;0]
plot(uv)
axis equal
%% With CHEBOP -- Lorenz
dom = [0 15];
N = chebop(@(t,u,v,w) [diff(u) - 10*(v - u);
    diff(v) - u.*(28 - w) + v;
    diff(w) - u.*v + (8/3)*w], dom);
N.lbc = @(u,v,w) [w - 20 ; v + 15; u + 14];
uvw = N\[0;0;0]
plot(uvw)
plot3(uvw{1},uvw{2}, uvw{3}, 'linewidth', 1.6), view(20,20)
axis([-20 20 -40 40 5 45]), grid on
xlabel 'x(t)', ylabel 'y(t)', zlabel 'z(t)'
title('A 3D Trajectory of the Lorenz Attractor - Chebfun solution', FS, 14)


