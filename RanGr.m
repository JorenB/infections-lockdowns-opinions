function [EdgeList] = RanGr(nNodes,p)
%RanGr generates a random network which has nNodes, and each edge exists
%with probability p, independent of existence of other edges
%   the function returns EdgeList which is an nNodesxn(Nodes)/1x3 table
%   where the first two columns contain nodes ids and the third entry is
%   assigned 0 or 1 depending on whether an edge between nodes exists
% the list EdgeList is ordered lexographically

% initialize the list

EdgeList=zeros(nNodes*(nNodes-1)/2,3);
 
% populate the list with node names
c=1;
for i1=1:nNodes-1
   for i2=(i1+1):nNodes
       EdgeList(c,1)=i1;
       EdgeList(c,2)=i2;
       c=c+1;
   end
end
%generate nNodes(nNodes-1)/2 uniform random numbers
rels=rand(nNodes*(nNodes-1)/2,1);
EdgeList(:,3)=rels<p;
end

