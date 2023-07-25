%This function uses matrix algebra to compute the unconditional
%variance-covariance matrix of a process X_{t+1}=A \tilde
%X_t+\epsilon, where \epsilon has variance-covariance matrix EN

% INPUTS: A: Matrix with eigenvalues less than one in absolute value
%         EN: positive-definite symmetric variance-covariance matrix 
%
% OUTPUTS: VarX0: The long-run variance matrix of the factors.


function VarX0 = Uncon_var(A, EN) 

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
