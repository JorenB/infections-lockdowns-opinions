function [Na,Nb]=NeighOp(People,FriendsIds)
%the function produces the number of A believers and B believers in the
%neighbourhoud of individual id
% inputs: 
% 1. table of people in the Population, People. The second column of
% the table contains id of each individual (1 - a, 0 - b).
% 2. table of direct ties of an individual (ids) in the
% opinion network
% output: number of direct ties with opinion a Na, number of direct ties with
% opinion b Nb

Na=sum(People(FriendsIds,2)==1);
Nb=numel(FriendsIds)-Na;
end
