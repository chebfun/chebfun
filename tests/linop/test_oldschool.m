function pass = test_oldschool

%% Test the 'old-school' linop syntaxes for ATAP compatibility
% Toby Driscoll, 3 March 2014

d = domain(-1,1);

%% Just test whether these will execute.
pass(1) = doesNotCrash( @() diff(d) );
pass(2) = doesNotCrash( @() eye(d) );
pass(3) = doesNotCrash( @() diag(chebfun('x',d([1 end])), d) );
pass(4) = doesNotCrash( @() cumsum(d) );
pass(5) = doesNotCrash( @() sum(d) );
pass(6) = doesNotCrash( @() feval(d, 0, 'left') );
pass(7) = doesNotCrash( @() zeros(d) );

%% Test the feval-style instantiation syntax.
D = diff(d);
pass(8) = norm( D(8) - matrix(D,8) ) < 2e-14;

%% Test boundary condition syntax
A = D^2;
A0 = matrix(A,10);   % version with no BCs
[z,e,s,dif] = linop.primitiveFunctionals(d([1 end]));
A = addConstraint(A,e(d(1)),0);
A = addConstraint(A,e(d(end))*D,0);
A1 = matrix(A,10);   % first two rows hold BCs
correct = A0; 
correct([1 end],:) = A1(1:2,:);  % classic row replacement
pass(9) = norm( correct - feval(A,10,'oldschool') ) < 2e-14;

end

function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME
    pass = false;
end
end
