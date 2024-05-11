% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function people = initInfect(people, initfInfected)
	% this  function takes a list of people People and the number of infected
	% individuals Ni and then randomly pick who will be infected, these are
	% picked will have their infection status (third entry of the respective
	% row) changed to numeric code denoting infectious status (1)
	% The function assumes that prior to its execution everybody has a
	% susceptible status (0)
	nPeople = size(people,1);
	nInfected = floor(initfInfected * nPeople);

	%generate the ids of these who hold opinion a
	ids = datasample(1:nPeople, nInfected, 'Replace', false);
	people(ids, 3) = 1;
end
