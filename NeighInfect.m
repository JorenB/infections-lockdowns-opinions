% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function infected=NeighInfect(People,FriendsIds)
%the function produces the number of infected in the passed list of direct
%ties in the physical network. The function assumes that the passed list is
%for a single person only

%extract infection status of direct ties
infected=sum(People(FriendsIds,3)==1);
end
