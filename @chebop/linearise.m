function [L, res, isLinear] = linearise(N, x, u, flag)
isLinear = true(1, 4);
numVars = nargin(N.op) - 1;

if ( nargin < 2 )
    x = chebfun(@(x) x, N.domain);
end
if ( nargin < 3 )
    % Initialise a zero ADCHEBFUN:
    zeroFun = chebfun(0, N.domain);
    u = cell(numVars, 1);
    for k = 1:numVars
        u{k} = zeroFun;
    end
    %                 u0 = chebmatrix(u0);
end

if ( isa(u, 'chebmatrix') )
    u = u.blocks;
end
if ( isa(u, 'chebfun') )
    u = {u};
end

for k = 1:numVars
    u{k} = seed(adchebfun(u{k}), k, numVars);
end


%%
% Evaluate the operators to get a linearisation:

% Experiment with using chebmatrices:
%             w = N.op(x, u{:});
%             if isa(w, 'adchebfun')
%                 L = w.jacobian;
%                 res = w.func;
%                 isLinear(1) = all(w.isConstant);
%             else
%                 L = cellfun(@(w) get(w, 'jacobian'), w.blocks, 'uniformoutput', false);
%                 L = vertcat(L{:});
%                 res = cellfun(@(w) get(w, 'func'), w.blocks, 'uniformoutput', false);
%                 res = chebmatrix(res);
%                 isLinear(1) = all(cellfun(@(w) get(w, 'isConstant'), w.blocks));
%             end

w = N.op(x, u{:});
L = linop(vertcat(get(w, 'jacobian')));
res = vertcat(get(w, 'func'));
isLinear(1) = all(vertcat(get(w, 'isConstant')));


% Merge the domains of L obtained from evaluating the operator part above, with
% the domain of N, as we want to respect the breakpoints originally assigned
% to N (when its domain was defined)
L.domain = chebfun.mergeDomains(L.domain, N.domain);

BC = linopConstraint();
%%
% Add BCs
if ( ~isempty(N.lbc) )
    lbcU = N.lbc(u{:});
    for k = 1:numel(lbcU)
        lbcUk = getElement(lbcU, k);
        lbcUk = feval(lbcUk, N.domain(1));
        if ( nargin == 4 ), lbcUk.func = -lbcUk.func; end
        BC = append(BC, lbcUk.jacobian, -lbcUk.func);
        %                     L = bc(L, lbcU(k).jacobian, -lbcU(k).func); %#ok<CPROP>
    end
    isLinear(2) = all(get(lbcU, 'isConstant'));
end

if ( ~isempty(N.rbc) )
    rbcU = N.rbc(u{:});
    for k = 1:numel(rbcU)
        rbcUk = getElement(rbcU, k);
        rbcUk = feval(rbcUk, N.domain(end));
        if ( nargin == 4 ), rbcUk.func = -rbcUk.func; end
        BC = append(BC, rbcUk.jacobian, -rbcUk.func);
        %                     L = bc(L, rbcU(k).jacobian, -rbcU(k).func); %#ok<CPROP>
    end
    isLinear(3) = all(get(rbcU, 'isConstant'));
end

if ( ~isempty(N.bc) )
    % Experimant with using chebmatrix to store adchebfuns:
    %                 bcU = N.bc(x, u{:});
    %                 con = cellfun(@(w) get(w, 'jacobian'), bcU.blocks, 'uniformoutput', false);
    %                 con = vertcat(con{:});
    %                 val = cellfun(@(w) get(w, 'func'), bcU.blocks, 'uniformoutput', false);
    %                 val = chebmatrix(val);
    %                 if ( nargin == 4 ), val = -val; end
    %                 L = bc(L, con, -val); %#ok<CPROP>
    
    bcU = N.bc(x, u{:});
    vals = cat(1, get(bcU, 'func'));
    if ( nargin == 4 ), vals = -vals; end
    for k = 1:numel(bcU)
        BC = append(BC, get(bcU, 'jacobian', k), -vals(k));
        %                     L = bc(L, bcU(k).jacobian, -vals(k)); %#ok<CPROP>
    end
    isLinear(4) = all(get(bcU, 'isConstant'));
end
L.constraint = BC;

end