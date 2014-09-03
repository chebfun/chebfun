function pass = test_Lorenz
%% With CHEBOP -- Lorenz
dom = [0 5];
N = chebop(@(t,u,v,w) [diff(u) - 10*(v - u);
    diff(v) - u.*(28 - w) + v;
    diff(w) - u.*v + (8/3)*w], dom);
N.lbc = @(u,v,w) [w - 20 ; v + 15; u + 14];
uvw = N\[0;0;0];
plot3(uvw{1},uvw{2}, uvw{3}, 'linewidth', 1.6), view(20,20)
axis([-20 20 -40 40 5 45]), grid on
xlabel 'x(t)', ylabel 'y(t)', zlabel 'z(t)'
title('A 3D Trajectory of the Lorenz Attractor - Chebfun solution', 'fontsize', 14)

pass = 1;
end