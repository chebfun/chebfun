function [u,lam,iter, retract] = newtonBVP(H,A,g,BCstruct,uinit, laminit, udot,lamdot)

% Find a tangent to the curve H(u,lambda)=0 at the given point, by solving
% a boundary-value problem.

bcFun = H.bc;
% if nargin < 4
%     bcFun = [];
% end


d = domain(udot);
x = chebfun(@(x) x,d);
% Begin by creating a chebarray for the constraint which goes at the
% bottom
Jm = [udot'*diag(sum(d)) lamdot];

% Convert lam to chebconst, and store udot as u for when we iterate
lam = laminit;   % to get autodiff data
% u = jacreset(uinit);
u = uinit;
% Compute Newton correction
accept = 0; iter = 0; retract = 0;
while ~accept
    % Evaluate residual, and linearise
    v = H.op(x, u, lam);             % should be zero, approximately!
    Au = A(u, lam);
    gu = g(u, lam);    
    
    % Create a chebarray
    S = linop([Au gu;Jm]);
    
    % rhs is the residual, with zero at the bottom for the functional
    % condition
    rhs = [-v;0];
    
    % Need to linearise and obtain superlinops to use for BCs if bcFun is not
    % empty
    if ~isempty(bcFun)
        % Evaluate bcFun
%         u = set(u,'funreturn',1);
        bcVals = bcFun(x,u,lam);
        % Differentiate each result to get a superlinop, consisting of a row
        % linop and a scalar
        
        % Assign values for the RHS of BCs
        S = addConstraint(S, BCstruct(1,:), 0);
        S = addConstraint(S, BCstruct(2,:), 0);
%         BCstruct.other(1).val = -0*bcVals(:,1).vals;
%         BCstruct.other(2).val = -0*bcVals(:,2).vals;
        
%         for bcCounter = 1:numel(bcVals)
%             rowlin = diff(bcVals(:,bcCounter),u,'linop');
%             rowsca = diff(bcVals(:,bcCounter),lam,'linop');
%             % Convert each part into a superlinop
%             rowop = [rowlin rowsca];
%             rowval = bcVals(:,bcCounter).vals;
%             BCstruct.other(bcCounter).op = rowop;
%             BCstruct.other(bcCounter).val = -rowval;
%         end
        % Revert to funreturn == 0
%         u = set(u,'funreturn',0);
    end
    
    
    dudlam = S\rhs;
    du = dudlam{1};
    dlam = dudlam{2};
    
    % Should be close to solution, so take full Newton steps
    u = u+du;
    lam = lam+dlam;
%     u = jacreset(u);
%     lam = jacreset(lam);
    
    
%     fprintf('Iter: %d, res. %6.4e .\n',iter, norm(du,2))
    iter = iter + 1;
    if norm(du,2) < 1e-5
        accept = 1;
    elseif iter >=10 % Too many iterations
        retract = 1; return
    end

end
end