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
A = [ f; g; h];
[normVal, col] = norm(A, 1);
pass(1) = (col == 3) && abs(normVal - (exp(1) - exp(-1))) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
                                   
pass(2) = abs(norm(A) - 2.372100421113536830) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
                                            
pass(3) = abs(norm(A, 'fro') - 2.372100421113536830) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
                                                                                      
[normVal, loc] = norm(A, inf);
pass(4) = (loc == 1) && abs(normVal - (exp(1) + sin(1) + cos(1))) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);

[normVal, loc] = norm(A, -inf);
pass(5) = (loc == -1) && abs(normVal - (exp(-1) + sin(1) + cos(1))) < ...
        (chebfun(A).vscale)*(chebfun(A).epslevel);
                                          
U = chebfun(@(x) [(1 + 0*x) exp(2*pi*1i*x) exp(2*pi*1i*2*x)], [0 1], pref);
S = diag([pi ; exp(1) ; 1]);
V = [1/sqrt(2) -1/sqrt(2) 0 ; 1/sqrt(2) 1/sqrt(2) 0 ; 0 0 1];
h = U*S*V';
A = [ h(:, 1); h(:, 2); h(:, 3) ];

pass(6) = abs(norm(A, 2) - pi) < (chebfun(A).vscale)...
                                 *(chebfun(A).epslevel);         
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