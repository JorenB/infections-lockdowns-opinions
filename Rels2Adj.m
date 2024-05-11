% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function AdjM=Rels2Adj(Rels,N)
% function that converts a reltionship matrix Rels to adjacency matrix in a
% population with N individuals, which are assumed to have sequential ids

AdjM=zeros(N,N);

num_rels=size(Rels,1);

for counter=1:num_rels
   id1=Rels(counter,1);
   id2=Rels(counter,2);
   AdjM(id1,id2)=1;
   AdjM(id2,id1)=1;
end