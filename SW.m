function A=SW(N,K,beta)
% This function constructs a  graph using Watts-Strogatz algorithm
%input: N: node cardinality, K: average degree, beta: rewiring
%probability parameter
%output: A: adjacency matrix of the graph 

ownName=mfilename;
Nodes=1:1:N;
if K<=log(N)
   error([ownName,': mean degree is insufficiently large']); 
elseif K>=N
    error([ownName,': mean degree is too large']);
end
if mod(K,2)~=0
    error([ownName,': K must be even integer']);
end
%adjacency matrix
A=zeros(N,N);
%create a lattice with fixed degree K
for id=1:N %id counter
    for count2=1:1:K/2
        %right
        if id+count2~=N
            A(id,mod(id+count2,N))=1;
            A(mod(id+count2,N),id)=1;    
        else
            A(id,N)=1;
            A(N,id)=1;
        end
        %left
        if id-count2~=0
            A(id,mod(id-count2,N))=1;
            A(mod(id-count2,N),id)=1;
        else
            A(id,N)=1;
            A(N,id)=1;
        end
        
    end
end

%initialize the probability of edge switching for K-2
for node=1:1:N
   switchprob=rand(1,K/2);
   switchN=switchprob<beta;
   %find partners
   target=find(A(node,:)~=0);
   % find these to whom its possible to reconnect
   NotPart=setdiff(Nodes,[node,Nodes(target)]);
   %locate all neighbors to avoid rewiring to existing partners
   %sort neighbors so the right most K/2 will be the first ones
   tmax=target(target>node);
   tmin=target(target<node);
   target=[tmax,tmin];
   target=target(1:K/2);
   
   for c=1:K/2
      t=target(c);
      %rewire the edge
      if switchN(c)>0
         %pull a random new target
         newtarget=unidrnd(numel(NotPart));
         %rewire
         %cancel previous partner
         A(node,t)=0;
         A(t,node)=0;
         %rewire to new partner
         A(node,NotPart(newtarget))=1;
         A(NotPart(newtarget),node)=1;
         %g=graph(A);
        % figure(31);plot(g);
        % pause(1);
         %update NotPart
         % these with whom node is connected now cannot be targeted
         NotPart=setdiff(NotPart,NotPart(newtarget));
         % these from whom node got disconnected can be connected to
         % again
         NotPart=sort([NotPart,t]);
      end
   end
end
end