% Copyright 2024 Alexandra Teslya & Joren Brunekreef
% This function calculates the average local clustering coefficient for a
% network given by the adjacency matrix A
function res=LCC(A)
    %check the input
    [m,N]=size(A);
    ownName=mfilename;
    if m~=N
        disp([ownName,': provided matrix is not square']);
    end
    if ~issymmetric(A)
        error([ownName,': Adjacency matrix is not symmetric']);
    end
    AAA=diag(A^3)/2;
    degree=sum(A,2);
    index=degree>1;
    LLCSW=zeros(N,1);
    LLCSW(index)=(2*AAA(index))./(degree(index).*(degree(index)-1));
    res=mean(LLCSW);
end
