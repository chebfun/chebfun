function [x,w,v] = equationPoints(disc)
%EQUATIONPOINTS Points at which collocation is enforced.

[x,w,v] = colloc.points(disc,1);

end