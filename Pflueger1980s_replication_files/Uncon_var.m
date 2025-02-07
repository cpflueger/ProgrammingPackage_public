function VarX0 = Uncon_var(A, EN) 
% Uncon_var is a function that solves the long-run variance matrix
% for linear state space models. It requires that the unconditional
% convariance matrix var(X_0) satisfies the equation A Var(X_0) A' + E(N) =
% Var(X_0). Here A is matrix with state variable persistence parameters, 
% see equation (156) for definition. And E(N) is matrix with mean state 
% variable volatility.
%
% INPUTS: A: Matrix with state variable persistence parameters in plant
%            equation, or model equation (156).
%         EN: Matrix with mean state variable volatility (E(N) in model 
%             equation (171)); see model equations (172) and (173) for 
%             expression defining E(N).
%
% OUTPUTS: VarX0: The long-run variance matrix of the factors.

% Check the size of matrix A. Matrix A should be a square matrix.
if size(A,1) ~= size (A,2)
    error('Error: The transition matrix must be a N*N square matrix.')
end

% The function uses the fact that Var(X_0) can be solved using 
% (I-A otimes A)\v, where v is a column vector whose elements are taken
% column-wise from EN. A otimes A is the Kronecker tensor product of matrix
% A and A itself.
n     = size(A,1);

% Compute the Kronecker product of A and A.
Astar = kron(A,A);

% Form a vector whose elements are taken from EN column-wisedly. Then
% compute the elements of VarX0. Finally reshape the results from a vector
% to a matrix.
VarX0 = reshape((eye(n^2)-Astar)\reshape(EN,n^2,1), n, n);
end
