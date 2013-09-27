function C = horzcat(varargin)

args = {};
for k = 1:numel(varargin)
    if ( isa(varargin{k}, 'chebmatrix') )
        args = [args, varargin{k}.blocks];
    else
        args = [args, varargin(k)];
    end
end

C = chebmatrix( args );

end