function [GProws,GPcols] = PartitionTreeMatrixDhillon(matrix)
% Recursively partition the rows and columns of the matrix using the
% bipartite model proposed by Dhillon in "Co-clustering documents and words
% using Bipartite Spectral Graph Partitioning"
%
% Input
%   matrix      a matrix where each column represents a sample
%
% Output
%   GProws      a GraphPart object ==> partitioning using rows as samples
%   GPcols      a GraphPart object ==> partitioning using cols as samples
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminary stuff

% constants
[rows,cols] = size(matrix);
N = max([rows,cols]);
jmax = max(3*floor(log2(N)),4);


%%% TRACKING VARIABLES FOR THE ROWS

% order index
sslac = ind_class(rows);
ind_rows = zeros(rows,1,sslac);
ind_rows = (1:rows)';

% regionstarts
rs_rows = zeros(N+1,jmax,sslac);
rs_rows(1,:) = 1;
rs_rows(2,1) = rows+1;


%%% TRACKING VARIABLES FOR THE COLUMNS

% order index
sslac = ind_class(cols);
ind_cols = zeros(cols,1,sslac);
ind_cols = (1:cols)';

% regionstarts
rs_cols = zeros(N+1,jmax,sslac);
rs_cols(1,:) = 1;
rs_cols(2,1) = cols+1;



%% 1. Partition the graph to yield rs, ind, and jmax

j = 1;
regioncount_rows = 0;
regioncount_cols = 0;

while regioncount_rows < rows || regioncount_cols < cols %regioncount < N
    % the number of true regions (i.e., has 1 or more nodes) on level j
    regioncount_rows = length(unique([rs_rows(:,j);0]))-2;
    regioncount_cols = length(unique([rs_cols(:,j);0]))-2;
    
    % the number of total regions (including those with no nodes)
    false_regioncount = max([nnz(rs_rows(:,j)), nnz(rs_cols(:,j))])-1;

    % add a column for level j+1, if necessary
    if j == jmax
        rs_rows(1,j+1) = 1;
        rs_cols(1,j+1) = 1;
        jmax = jmax+1;
    end

    % for tracking the child regions (for both rows and columns)
    rr = 1;

    % cycle through the parent regions
    for r = 1:false_regioncount
        % the start of the parent region and 1 node after the end of
        % the parent region
        rs1_rows = rs_rows(r,j);
        rs2_rows = rs_rows(r+1,j);
        rs1_cols = rs_cols(r,j);
        rs2_cols = rs_cols(r+1,j);

        % the number of nodes in the parent region
        n_rows = double(rs2_rows-rs1_rows);
        n_cols = double(rs2_cols-rs1_cols);
        
        if n_rows <= 1 && n_cols > 1
            indrs_rows = ind_rows(rs1_rows);
            indrs_cols = ind_cols(rs1_cols:rs2_cols-1);
            
            % compute the left and right singular vectors
            [~,v_cols] = second_largest_singular_vectors(matrix(indrs_rows,indrs_cols));

            % partition the current region
            pm_cols = PartitionFiedler(1,1,v_cols);
            
            % determine the number of points in child region 1
            n1_cols = sum(pm_cols > 0);
            
            % update the column indexing
            ind_cols(rs1_cols:rs1_cols+n1_cols-1) = indrs_cols(pm_cols > 0);
            ind_cols(rs1_cols+n1_cols:rs2_cols-1) = indrs_cols(pm_cols < 0);

            % update the row region tracking
            rs_rows(rr+1,j+1) = rs1_rows;
            rs_rows(rr+2,j+1) = rs2_rows;
            
            % update the column region tracking
            rs_cols(rr+1,j+1) = rs1_cols+n1_cols;
            rs_cols(rr+2,j+1) = rs2_cols;
            
            rr = rr+2;
            
        elseif n_rows > 1 && n_cols <= 1
            indrs_rows = ind_rows(rs1_rows:rs2_rows-1);
            indrs_cols = ind_cols(rs1_cols);
            
            % compute the left and right singular vectors
            [u_rows,~] = second_largest_singular_vectors(matrix(indrs_rows,indrs_cols));

            % partition the current region
            pm_rows = PartitionFiedler(1,1,u_rows);
            
            % determine the number of points in child region 1
            n1_rows = sum(pm_rows > 0);     
            
            % update the row indexing
            ind_rows(rs1_rows:rs1_rows+n1_rows-1) = indrs_rows(pm_rows > 0);
            ind_rows(rs1_rows+n1_rows:rs2_rows-1) = indrs_rows(pm_rows < 0);

            % update the row region tracking
            rs_rows(rr+1,j+1) = rs1_rows+n1_rows;
            rs_rows(rr+2,j+1) = rs2_rows;
            
            % update the column region tracking
            rs_cols(rr+1,j+1) = rs1_cols;
            rs_cols(rr+2,j+1) = rs2_cols;
            
            rr = rr+2;
            
        elseif n_rows > 1 && n_cols > 1
            indrs_rows = ind_rows(rs1_rows:rs2_rows-1);
            indrs_cols = ind_cols(rs1_cols:rs2_cols-1);
            
            % compute the left and right singular vectors
            [u_rows,v_cols] = second_largest_singular_vectors(matrix(indrs_rows,indrs_cols));

            % partition the current region
            pm = PartitionFiedler(1,1,[u_rows; v_cols]);
            pm_rows = pm(1:n_rows);
            pm_cols = pm(n_rows+1:end);
            
            % determine the number of points in child region 1
            n1_rows = sum(pm_rows > 0);
            n1_cols = sum(pm_cols > 0);
            
            % update the row indexing
            ind_rows(rs1_rows:rs1_rows+n1_rows-1) = indrs_rows(pm_rows > 0);
            ind_rows(rs1_rows+n1_rows:rs2_rows-1) = indrs_rows(pm_rows < 0);

            % update the column indexing
            ind_cols(rs1_cols:rs1_cols+n1_cols-1) = indrs_cols(pm_cols > 0);
            ind_cols(rs1_cols+n1_cols:rs2_cols-1) = indrs_cols(pm_cols < 0);

            % update the row region tracking
            rs_rows(rr+1,j+1) = rs1_rows+n1_rows;
            rs_rows(rr+2,j+1) = rs2_rows;
            
            % update the column region tracking
            rs_cols(rr+1,j+1) = rs1_cols+n1_cols;
            rs_cols(rr+2,j+1) = rs2_cols;
            
            rr = rr+2;
            
        else
            rs_rows(rr+1,j+1) = rs2_rows;
            rs_cols(rr+1,j+1) = rs2_cols;
            rr = rr+1;
        end
    end

    j = j+1;
end


% get rid of excess columns in rs_rows and rs_cols
rs_rows(:,j:end) = [];
rs_cols(:,j:end) = [];
[~,jmax] = size(rs_rows);

% remove duplicates
[M_rows,~] = size(rs_rows);
[M_cols,~] = size(rs_cols);
for j = jmax:-1:1
    %%% ROWS
    % remove duplicate regionstarts
    unique_rows = unique(rs_rows(:,j));
    if unique_rows(1) == 0
        unique_rows(1) = [];
    end
    rs_rows(:,j) = [unique_rows; zeros(M_rows-length(unique_rows),1)];
    
    % remove duplicate columns of rs
    if j < jmax && sum(rs_rows(:,j) ~= rs_rows(:,j+1)) == 0
        rs_rows(:,j+1) = [];
    end
        
    
    %%% COLUMNS
    % remove duplicate regionstarts
    unique_cols = unique(rs_cols(:,j));
    if unique_cols(1) == 0
        unique_cols(1) = [];
    end
    rs_cols(:,j) = [unique_cols; zeros(M_cols-length(unique_cols),1)];
    
    % remove duplicate columns of rs
    if j < jmax && sum(rs_cols(:,j) ~= rs_cols(:,j+1)) == 0
        rs_cols(:,j+1) = [];
    end
end

% remove unnecessary rows and columns in the regionstarts
if rows < M_rows-1
    rs_rows(rows+2:end,:) = [];
end
if cols < M_cols-1
    rs_cols(cols+2:end,:) = [];
end

% create GraphPart objects
GProws = GraphPart(ind_rows,rs_rows,[],[],[],[],[],'Dhillon');
GPcols = GraphPart(ind_cols,rs_cols,[],[],[],[],[],'Dhillon');

% make sure the GraphPart objects meet our requirements
GProws = GP_Fix(GProws,GraphSig(sparse(rows,rows)));
GPcols = GP_Fix(GPcols,GraphSig(sparse(cols,cols)));


end






function [u,v] = second_largest_singular_vectors(A)
% Compute the second largest left and right singular vectors of 
% An = D1^-0.5 A D2^-0.5

[rows,cols] = size(A);

% compute D1 and D2 and make sure there are no zeros
D1 = sum(A,2);
D1(D1 == 0) = max([0.01, min(D1(D1 > 0))/10]);
D2 = sum(A,1);
D2(D2 == 0) = max([0.01, min(D2(D2 > 0))/10]);

% compute the singular vectors
try
    if rows <= 128 || cols <= 128
        [u,~,v] = svd( full( bsxfun(@times, bsxfun(@times, (D1).^(-0.5), A), (D2).^(-0.5) ) ) );
    else
        [u,~,v] = svds(bsxfun(@times, bsxfun(@times, (D1).^(-0.5), A), (D2).^(-0.5) ), 2);
    end
    u = u(:,2);
    v = v(:,2);
catch
    [u,~,v] = svd( full( bsxfun(@times, bsxfun(@times, (D1).^(-0.5), A), (D2).^(-0.5) ) ) );
    
    % extract the second left singular vector
    [~,ucols] = size(u);
    if ucols == 1
        u = u(:,end);
    else
        u = u(:,end-1);
    end
    
    % extract the second right singular vector
    [~,vcols] = size(v);
    if vcols == 1
        v = v(:,end);
    else
        v = v(:,end-1);
    end
end

% extract the 2nd singular vectors and multiply by D_i^-0.5
u = D1.^-0.5 .* u;
v = (D2').^-0.5 .* v;
end