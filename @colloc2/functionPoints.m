function [x,w] = functionPoints(disc)
%FUNCTIONPOINTS Points at which functions are discretized.

% In COLLOC2, functions are discretized at 2nd kind points but equations
% are enforced at 1st kind points, to avoid duplication at boundaries.

[x,w] = colloc.points(disc,2);

end