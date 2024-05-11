% Copyright 2024 Alexandra Teslya & Joren Brunekreef
% This function calculates the average path length for a
% network given by the adjacency matrix A
function res=APL(A)
    %check the input
    [m,N]=size(A);
    ownName=mfilename;
    if m~=N
        disp([ownName,': Adjacency matrix is not square']);
    end
    if ~issymmetric(A)
        error([ownName,': Adjacency matrix is not symmetric']);
    end
    G=graph(A);
    d=distances(G);
    res=sum(sum(d))/(N*(N-1));
end