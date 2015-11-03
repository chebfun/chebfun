function pass = test_cellOperator
% TEST_CELLOPERATOR   Test a chebop whose OP is a function that returns a cell

N = chebop(@(r,u) odefun(r,u), [-1 1]);
N.lbc = @(u) [u{1} - 1 ; u{2}];
N.rbc = @(u) [u{1} ; u{2} - 1];
% Need to specify the number of input variables, as otherwise, we won't know the
% number of variables involved:
N.numVars = 2; 
u = N\0;
% Did we get the results expected?
pass = norm(N(u)) < 1e-10 && norm(feval(N.lbc(u), -1)) < 1e-13 && ...
    norm(feval(N.rbc(u), 1)) < 1e-13;

    function op = odefun(r, u)
        for n = 1:2
            op{n} = diff(u{n}, 2) + n*u{n};
        end
    end

end