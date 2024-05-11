% Copyright 2024 Alexandra Teslya & Joren Brunekreef
% convert adjacency matrix to relationship list of the form [id1,id2]
% the graph inducing the adjecency matrix does not have node self-loops
% Theresulting list Rels is ordered lexographically
function Rels = Adj2Rels(AdjM)
    % adjacency matrices are assumed to be in sparse format
    A=full(AdjM);
    %check input
    [N,M]=size(A);
    ownName=mfilename;
    if N~=M
       error([ownName,': Matrix AdjM has to be a square matrix']);
    end

    if ~issymmetric(AdjM)
        error([ownName,': Matrix AdjM has to be symmetric']);
    end
    
    % initialize the overall contained, could be more efficient, as the
    % size of Rels is sum(sum(AdjM))/2 x 2
    Rels=[];
    
    for node1=1:N-1
        nodes=node1+1:1:N;
        ind=A(node1,node1+1:N)>0;
        nrels=sum(ind); % number of relationships of node1
        rels(:,1)=zeros(1,nrels)+node1;
        rels(:,2)=nodes(ind);
        Rels=[Rels;rels];
        clear rels;
    end
end