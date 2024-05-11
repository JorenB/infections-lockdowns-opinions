% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function [t, y, people] = ssaBurn(people, opNetRels, physNetRels, cfg)
	% this function generates burn in dynamics for the opinion dynamics and
	% stops when the time counter reaches the value of cfg.tBurn
	% input: table of individuals people, opinion network OpNet, physical
	% network PhysNet, parameter list pars, time of the simulation tBurn,
	% number of points where measurements of the networks will be taken

	%define threshold number of steps
	nsteps=10000;
	nPeople = cfg.nPeople;

	%people=[id,op,ep,Na,Nd,Ninfect];id=1..N,op=0,1,ep=0,1,2


	%initialize y=[sa,ia,ra,ida,sb,ib,rb,idb]

	peopleA = people(find(people(:,2) == 1),:);
	peopleB = people(find(people(:,2) == 0),:);
	Sa = peopleA(find(peopleA(:,3) == 0),:);
	Ia = peopleA(find(peopleA(:,3) == 1),:);
	Ra = peopleA(find(peopleA(:,3) == 2),:);
	Sb = peopleB(find(peopleB(:,3) == 0),:);
	Ib = peopleB(find(peopleB(:,3) == 1),:);
	Rb = peopleB(find(peopleB(:,3) == 2),:);
	y = cellfun(@(x) size(x, 1), {Sa Ia Ra Sb Ib Rb});

	% initialize time array
	t = [];
	t(1) = 0;


	% check that there are no infected people before the opinion burn
	assert(size(Ia, 1) + size(Ib, 1) == 0);

	counter = 1;

	Adj = Rels2Adj(opNetRels, nPeople);

	%create friends lists: opinion network
	for id = 1:nPeople
		FriendsOp{id} = find(Adj(id,:) > 0);

		% for each individual collect their Na and Nb neighbours
		[Na, Nb] = NeighOp(people, FriendsOp{id}); 
		people(id, 4) = Na;
		people(id, 5) = Nb;
	end

	%calculate opinion switch propensities with nInfected = 0 and ld_duration = 0
	osps = opSwitchPropensities(people, 1:nPeople, cfg, 0, 0);

	osps = [osps, zeros(1, nPeople)];

	sumOsps = sum(osps);
	if sumOsps > 0
		fl = 0;
	else
		fl = 1;
	end

	tic;
	while (t(counter) < cfg.tBurn) & ~fl
		%disp('b1') % begin burn, determine event type/id, update people and network
		rand1 = rand;
		rand2 = rand;

		%determine the time to the event
		deltaT = log(1 / rand1) / sumOsps;

		t(counter + 1) = t(counter) + deltaT;
		cumSumOsps = cumsum(osps);
		j = find(rand2 * sumOsps < cumSumOsps, 1);

		id = j;
		opOld = people(id, 2);
		opNew = 1 - opOld;
		infected = people(id, 3);
		people(id, 2) = opNew;


		IdsOp = [id];
		idNeigh = FriendsOp{id}; %unique(sort([OpNet(find(OpNet(:,1) == id),2);OpNet(find(OpNet(:,2) == id),1)]));
		%update status of id, update Na and Nb of its neighbours

		for idN = idNeigh
			[Na,Nb] = NeighOp(people, FriendsOp{idN});
			people(idN, 4) = Na;
			people(idN, 5) = Nb;
		end

		IdsOp = [IdsOp idNeigh];
		%record in y
		y(counter + 1,:) = y(counter,:);
		y(counter + 1, infected + 1 + 3 * (1 - opNew)) = y(counter + 1, infected + 1 + 3 * (1-opNew)) + 1;
		y(counter + 1, infected + 1 + 3 * opNew) = y(counter + 1, infected + 1 + 3 * opNew) - 1;

		osps(IdsOp) = opSwitchPropensities(people, IdsOp, cfg, 0, 0);

		counter = counter + 1;
		sumOsps = sum(osps);
		if sumOsps == 0
			fl = 1;
		end
	end
	toc

	%check whether we got to the final time, and if not 'bookend' the
	%trajectory'
	if t(counter) < cfg.tBurn
		t(counter + 1) = cfg.tBurn;
		y(counter + 1, :) = y(counter, :);
	end
end
