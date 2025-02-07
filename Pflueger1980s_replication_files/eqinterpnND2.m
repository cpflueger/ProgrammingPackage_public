%% Function to interpolate asset prices in simulation 
%griddedInterpolant is the N-dimensional interpolation routine available in
%MATLAB
function Vq = eqinterpnND2(varargin)
X         = varargin(1:5);
V         = varargin{6};
Xq        = varargin(7:11);
CTE       = varargin(12);
for i=1:length(V)
    F     = griddedInterpolant(X, V{i}, 'linear','none');
    Vq(i) = exp(F(Xq{:})+CTE{1});
end
return
end
