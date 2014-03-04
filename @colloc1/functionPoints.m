function [x,w,v] = functionPoints(disc)
%FUNCTIONPOINTS Points at which functions are discretized.

% In COLLOC1, functions are discretized at 1nd kind points and equations
% are enforced at 1st kind points.

[x,w,v] = colloc.points(disc,1);

end