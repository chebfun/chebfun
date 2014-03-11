function [x,w,v] = functionPoints(disc)
%FUNCTIONPOINTS Points at which functions are discretized.

[x,w,v] = colloc.points(disc,1);

end