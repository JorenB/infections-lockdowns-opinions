% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function num=ClusteringCoef(A)
% this function takes as the input an adjacency matrix and outputs the 
% global clustering coefficient equal to 3
% x number of closed triangles/number of triplets

nameStr = mfilename;

if size(A,1)~=size(A,2)
    error([nameStr,': Adjacency matrix is not square, it is ',num2str(size(A,1)),'x',num2str(size(A,2))]);
end

if ~issymmetric(A)
    error([nameStr,': Adjacency matrix is not symmetric']);
end
ntriangl=trace(A^3);
A2=A^2;
ntriples=sum(sum(A2))-trace(A2);
if ntriangl==0
    num=0;
else
    num=ntriangl/ntriples;
end
end