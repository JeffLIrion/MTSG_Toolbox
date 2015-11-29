function [GProws,GPcols,Wrows,Wcols] = PartitionTreeMatrixFiedler(matrix,weightfun,iter)
% Generate partition trees for the rows and columns of a matrix by forming 
% weight matrices and partitioning them via Fiedler vectors.  
%
% Input
%   matrix      a matrix where each column represents a sample
%   weightfun   the function for computing the weight between two nodes (if
%               not specified, inverse Euclidean weights are used)
%   iter        the number of iterations to perform
%
% Output
%   GProws      a GraphPart object ==> partitioning using rows as samples
%   GPcols      a GraphPart object ==> partitioning using cols as samples
%   Wrows       the weight matrix on the rows
%   Wcols       the weight matrix on the columns
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
% Wed Oct 14 21:06:07 2015  added by  Naoki Saito



% the number of iterations used
if ~exist('iter','var')
    iter = 10;
end

% the weight function to be used
if ~exist('weightfun','var')
    [weightfun,use_sqrt,use_RegInvEuc] = weight_function_matrix([]);
else
    [weightfun,use_sqrt,use_RegInvEuc] = weight_function_matrix(weightfun);
end

% define the non-diagonal indices if needed
if use_sqrt || use_RegInvEuc
    [rows,cols] = size(matrix);
    non_diags_rows = setdiff(1:rows^2, 1:rows+1:rows^2);
    non_diags_cols = setdiff(1:cols^2, 1:cols+1:cols^2);
end

% repeatedly partition the columns and rows
for i = 1:iter
    %%% Construct a partition tree on the columns of the matrix
    if i == 1
        Wcols = DualWeightMatrix(matrix,1,weightfun);
    else
        Wcols = DualWeightMatrix(matrix,GProws,weightfun);
    end
    
    % modify the weight matrix
    if use_sqrt
        Wcols = modify_sqrt(Wcols,non_diags_cols);
    elseif use_RegInvEuc
        Wcols = modify_RegInvEuc(Wcols,non_diags_cols);
    end
    
    % generate the recursive partitioning
    GPcols = PartitionTreeFiedlerSweep(GraphSig(Wcols),1);
%     GPcols = PartitionTreeFiedler(GraphSig(Wcols),1);
%     GPcols = PartitionTreeFiedlerkmeans(Gcols,1);
%     GPcols = PartitionTreeFiedlerMedian(Gcols,1);
%     GPcols = PartitionTreeFiedlerSweep(Gcols,1);


    %%% Construct a partition tree on the rows of the matrix
    if i == 1
        Wrows = DualWeightMatrix(matrix',1,weightfun);
    else
        Wrows = DualWeightMatrix(matrix',GPcols,weightfun);
    end
    
    % modify the weight matrix
    if use_sqrt
        Wrows = modify_sqrt(Wrows,non_diags_rows);
    elseif use_RegInvEuc
        Wrows = modify_RegInvEuc(Wrows,non_diags_rows);
    end
    
    % generate the recursive partitioning
    GProws = PartitionTreeFiedlerSweep(GraphSig(Wrows),1);
%     GProws = PartitionTreeFiedler(GraphSig(Wrows),1);
%     GProws = PartitionTreeFiedlerkmeans(Grows,1);
%     GProws = PartitionTreeFiedlerMedian(Grows,1);
%     GProws = PartitionTreeFiedlerSweep(Grows,1);
    

end


% generate the final weight matrix for the columns, which is based on the
% final row partitioning
Wcols = DualWeightMatrix(matrix,GProws,weightfun);

% modify the weight matrix
if use_sqrt
    Wcols = modify_sqrt(Wcols,non_diags_cols);
elseif use_RegInvEuc
    Wcols = modify_RegInvEuc(Wcols,non_diags_cols);
end


end



function W = modify_sqrt(W,non_diags)

N = length(W);
Wmin = min(W(non_diags));
Wmax = max(W(non_diags));
W = 2*(Wmax-W)/(Wmax-Wmin)-1; % linear
W = 0.5 + 0.5*sign(W).*nthroot(abs(W),2); % s-shape
W(1:1+N:N^2) = 0;

end


function W = modify_RegInvEuc(W,non_diags)

N = length(W);
Wmin = min(W(non_diags));
Wmax = max(W(non_diags));
W(1:1+N:N^2) = Wmin;
if isinf(Wmax)
    Wmax = 10*max(W(~isinf(W)))-9*Wmin;
    W(isinf(W)) = Wmax;
end

if abs(Wmin-Wmax) < 10^-13
    W = ones(N)-eye(N);
else
    W = (W-Wmin)/(Wmax-Wmin);
end

end