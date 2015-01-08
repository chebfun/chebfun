function [t,tau] = tangentBVP(H,A,g,BCstruct,u,lambda,told,tauold)

% Find a tangent to the curve H(u,lambda)=0 at the given point, by solving
% a boundary-value problem.

bcFun = H.bc;

% if nargin < 7
%     bcFun = [];
% end

d = domain(u);
x = chebfun(@(x) x, d);

% Differentiate to get H' = [A | f ]. 
lam = lambda;   % to get autodiff data
v = H.op(x,u,lam);              % should be zero, approximately!

Au = A(u,lam);
gu = g(u,lam);
% A = diff(v,u,'linop');
% g = diff(v,lam);

% u = set(u,'funreturn',1);

% Create a superlinop for the constraint which goes at the bottom
Jm = [told'*diag(sum(d)) tauold];

% Jul = J(u,lam);

% Jlin = diff(Jul,u,'linop');
% m = 1;%diff(Jul,lam);
% m = diff(Jul,lam);
S = linop([Au gu;Jm]);

% Need to linearise and obtain superlinops to use for BCs if bcFun is not
% empty
if ~isempty(bcFun)
    % Evaluate bcFun
%     bcVals = bcFun(x,u,lam);
    % Differentiate each result to get a superlinop, consisting of a row
    % linop and a scalar
        
    % Assign rhs of BCs. The op fields have been passed down.
    S = addConstraint(S, BCstruct(1,:), 0);
    S = addConstraint(S, BCstruct(2,:), 0);
    
%     for bcCounter = 1:numel(bcVals)
%        rowlin = diff(bcVals(:,bcCounter),u,'linop');
%        rowsca = diff(bcVals(:,bcCounter),lam,'linop');
%        % Convert each part into a superlinop
%        rowop = [rowlin rowsca];
%        rowval = bcVals(:,bcCounter).vals;
%        BCstruct.other(bcCounter).op = rowop;
%        BCstruct.other(bcCounter).val = 0; % Tangent problem always has conditions == 0
%     end
end
% BCstruct(2).val = -1;
rhs = [chebfun(0,d);1];

ttau = S\rhs;

% Extract function and scalar
t = ttau{1};
tau = ttau{2};

scale = sqrt(t'*t + tau^2);

t = t/scale;
tau = tau/scale;

end