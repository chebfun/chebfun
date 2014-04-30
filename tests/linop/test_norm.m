function pass = test_norm(pref)
% HM, 30 Apr 2014

%% 
% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

%% 
% The entries of A are only CHEBFUN or DOUBLE.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) cos(x), [-1 -0.5 0 0.5 1], pref);
h = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
A = [ f; g; h;];

pass(1) = abs(norm(A) - 2.372100421113536830) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
                                            
pass(2) = abs(norm(A, 'fro') - 2.372100421113536830) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
    
pass(3) = abs(norm(A, 2) - 2.372100421113536830) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
    
%% [TODO]
% A has entries of all types: OPERATORBLOCK, FUNCTIONBLOCK,
% CHEBFUN and DOUBLE. 
% d = [-2 2];                   % function domain
% I = operatorBlock.eye(d);     % identity operator
% D = operatorBlock.diff(d);    % differentiation operator
% x = chebfun(@(x) x, d);       % the variable x on d
% M = operatorBlock.mult(x.^2); % multiplication operator
% S = functionalBlock.sum(d);   % integration functional 
% E = functionalBlock.eval(d);  % evaluation functional generator
% A = [ I+D, abs(x), M; S, 0, E(2); D, x.^2, I ];
% [normOfA, loc] = norm(A);

end