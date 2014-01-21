function [order1, order2] = taylorTesting(f,hMax,plotting)

% By default, plotting is off.
if ( nargin < 3 )
    plotting = 0;
end

% Seed random generator to ensure same values.
seedRNG(6179);

% Initial u we evaluate f at
u = adchebfun(0.1*rand(8,1)+.5,[-1 1]);

% Chebfun used to create a perturbation
p = chebfun(0.01*rand(8,1)+.05,[-1 1]);

%% Taylor testing computations

% Vectors for storing information about convergence
nDiff1 = zeros(hMax,1);
nDiff2 = zeros(hMax,1);

% Factor we decrease by
fact = 5;

% Evaluate f at u
fu = f(u);

% Extract its Frechet derivative
fuJac = fu.jacobian;

for hCounter = 1:hMax
    % Compute a perturbation function
    pert = fact^-(hCounter+1)*p;

    % Evaluate f at the function obtained by adding the perturbation to u
    fuPert = f(u+pert);
    
    % Compute information for Taylor testing.
    
    % Expect this to be O(h)
    nDiff1(hCounter) = norm(fuPert - fu);
    
    % Use derivative information as well. Expect this to be O(h^2)
    nDiff2(hCounter) = norm(fuPert - fu - fuJac*pert);
end

%% Compute the orders of convergence
% Expect this to have entries close to 1
order1 = diff(log(nDiff1))./diff(log(fact.^-(2:hMax+1))');
% Expect this to have entries close to 2
order2 = diff(log(nDiff2))./diff(log(fact.^-(2:hMax+1))');


%% Plotting

if ( plotting == 1 )
    loglog(fact.^-(1:hMax),nDiff1,'-*'), hold on
    % Expected rate of convergence
    loglog([1 5^-20],[1 5^-20],'k--')
    
    loglog(fact.^-(1:hMax),nDiff2,'r-*'), hold on
    % Expected rate of convergence
    loglog([1 5^-10],[1 5^-20],'k--')
    set(gca,'XDir','reverse'), shg
    grid on
    shg
    axis equal
end

end