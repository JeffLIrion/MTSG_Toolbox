function G = Ggrid(Nx,Ny,~,~)
% Generate a GraphSig object for a 2-D grid that is Nx by Ny.  
%
% Input
%   Nx      the number of points in the x direction
%   Ny      the number of points in the y direction
%   ~       if a 3rd argument is given, use diagonal edges as well
%   ~       if a 4th argument is given, use the complete graph
%
% Output
%   G       the GraphSig ojbect
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the total number of nodes
N = Nx*Ny;

% the xy coordinates
xy = [(1:N)', repmat((1:Ny)',Nx,1)];
xy(:,1) = ceil(xy(:,1)/Ny);


if nargin == 2
    % make 1-D weight matrices
    ex = ones(Nx,1);
    Wx = spdiags([ ex 0*ex ex], [-1 0 1],Nx,Nx);
    ey = ones(Ny,1);
    Wy = spdiags([ ey 0*ey ey], [-1 0 1],Ny,Ny);

    % form the 2-D weight matrix
    W = full(kron(speye(Nx),Wy) + kron(Wx,speye(Ny)));   
else
    W = MakeDistMatrix(xy);
    W(1:N+1:N*N) = 10;
    W = W.^-1;
    
    if nargin == 3
        % keep only left, right, up, down, and diagonal edges
        W(W < 0.7) = 0;
    else
        % set the diagonal to zero (i.e., no loops)
        W(1:N+1:N*N) = 0;
    end
end


f = (1:N)';
G = GraphSig(W,xy,f,sprintf('%d by %d grid',Nx,Ny));


end
