function pass = test_jump_scaled
% TEST_JUMP_SCALED   Test a scaled jump. See #1694

A = chebop(@(x,u) diff(u,2) - u+x, [-1,1]);
A.lbc = .2;
A.rbc = 0;
A.bc = @(x,u) [ u(.1, 'left') - 4*u(.1, 'right') - 2.2; % Scaled jump in sol
    jump(diff(u),.1,1)];    % Continuous first derivative
u = A\0;
pass = abs(u(.1,'left') - 4*u(.1,'right') - 2.2 ) < 1e-10;

end