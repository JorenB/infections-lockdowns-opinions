function [tRec, yRec, ldsRec, nTransmissions] = ssaHybrid(people, opNetRels, physNetRels, cfg)
	%People=[id,op,ep,Na,Nd,Ninfect];id=1..N,op=0,1,ep=0,1,2

	peopleA = people(find(people(:,2) == 1),:);
	peopleB = people(find(people(:,2) == 0),:);
	Sa = peopleA(find(peopleA(:,3) == 0),:);
	Ia = peopleA(find(peopleA(:,3) == 1),:);
	Ra = peopleA(find(peopleA(:,3) == 2),:);
	Sb = peopleB(find(peopleB(:,3) == 0),:);
	Ib = peopleB(find(peopleB(:,3) == 1),:);
	Rb = peopleB(find(peopleB(:,3) == 2),:);
	y = cellfun(@(x) size(x, 1), {Sa Ia Ra Sb Ib Rb});


	yRec = cellfun(@(x) size(x, 1), {Sa Ia Ra Sb Ib Rb});

	nTransmissions = zeros(1, cfg.nPeople);
	nTransmissions(find(people(:, 3) == 1)) = 1;

	nInfected = y(:, 2) + y(:, 5);

	% start without lockdown
	ld = 0;
	ld_duration = 0;

	% initialize time array
	t = [];
	t(1) = 0;

	lds = [0];
	ldsRec = [0];
	
	counter = 1;
	tRes = 0.01;
	tCur = tRes;
	tRec = [0];

	tic
	AdjOp = Rels2Adj(opNetRels, cfg.nPeople);

	AdjPhys = Rels2Adj(physNetRels, cfg.nPeople);

	toc

	
	fprintf("indexing friends list\n");
	tic;
	FriendsOp = {};
	FriendsPhys = {};
	%create friends lists: opinion network & Physical
	for id = 1:cfg.nPeople
		FriendsOp{id} = find(AdjOp(id,:) > 0);
		FriendsPhys{id} = find(AdjPhys(id, :) > 0);

		% for each individual collect numbers of opinion/disease status neighbours
		[nOp_a, nOp_b] = NeighOp(people, FriendsOp{id}); 
		nInfectedNeighbors = NeighInfect(people, FriendsPhys{id}); 
		people(id, 4) = nOp_a;
		people(id, 5) = nOp_b;
		people(id, 6) = nInfectedNeighbors;
		people(id, 7) = numel(FriendsPhys{id}) - nInfectedNeighbors;
	end
	toc

	fprintf("computing propensities\n");
	osps = opSwitchPropensities(people, 1:cfg.nPeople, cfg, nInfected, ld_duration);
	%sprintf('%0.6f ', osps)

	% if more than zero people are infected
	if y(counter, 2) + y(counter, 5) > 0
		osps(cfg.nPeople + 1 : 2 * cfg.nPeople) = infectionPropensities(people, opNetRels, physNetRels, 1:cfg.nPeople, cfg);
	else
		osps = [osps, zeros(1, cfg.nPeople)];
	end

	sumOsps = sum(osps);
	if sumOsps > 0
		fl = 0;
	else
		fl = 1;
	end

	toc

	tic
	% main ssa loop
	fprintf("ssa start\n");
	tfc = 0;
	while (t(counter) < cfg.T_f) & ~fl
		rand1 = rand;
		rand2 = rand;

		sumOsps = sum(osps);

		%determine the time to the event
		deltaT = log(1 / rand1) / sumOsps;
		assert(deltaT >= 0); % seems superfluous

		t(counter + 1) = t(counter) + deltaT;
		cumSumOsps = cumsum(osps);
		j = find(rand2 * sumOsps < cumSumOsps, 1);
		if j <= cfg.nPeople % opinion switch
			id = j;
			opOld = people(id, 2);
			opNew = 1 - opOld;
			infStatus = people(id, 3);
			people(id, 2) = opNew;

			% update Adjacency matrices for Na and Nb
			% coefficients
			if opOld == 1 % old opinion was a
				
				%fprintf("A -> B\n");
				if ld
					inds1 = find(physNetRels(:, 1) == id);
					inds2 = find(physNetRels(:, 2) == id);
					inds = unique([inds1; inds2])';

					reinstates = 0;
					for ind = inds
						id1 = physNetRels(ind, 1);
						id2 = physNetRels(ind, 2);
						reinstateProb = 0;
						% both peers are now of opinion B
						if people(id1, 2) + people(id2, 2) == 0
							reinstateProb = 1 - (1 - cfg.alpha) * cfg.q;
						else
							reinstateProb = 1 - cfg.q;
						end

						if rand < reinstateProb
							physNetRels(ind, 3) = cfg.c_phys;
							reinstates = reinstates + 1;
						end
					end
					%fprintf("reinstates: %d / %d\n", reinstates, size(inds, 2));
				end

			else % old opinion was b
				
				%fprintf("B -> A\n");
				if ld
					inds1 = find(physNetRels(:, 1) == id);
					inds2 = find(physNetRels(:, 2) == id);
					inds = unique([inds1; inds2])';

					cuts = 0;
					for ind = inds
						id1 = physNetRels(ind, 1);
						id2 = physNetRels(ind, 2);

						% at least one peer is necessarily of opinion A
						cutProb = cfg.q;

						% we do not reinstate edges when switching from B to A
						if rand < cutProb
							physNetRels(ind, 3) = 0;
							cuts = cuts + 1;
						end
					end
					%fprintf("cuts: %d / %d\n", cuts, size(inds, 2));
				end
			end

			%collect ids of everybody who needs to be updated
			IdsOp = [id];
			idNeigh = FriendsOp{id};
			%update status of id, update Na and Nb of its neighbours

			for idN = idNeigh
				[Na,Nb] = NeighOp(people, FriendsOp{idN});
				people(idN, 4) = Na;
				people(idN, 5) = Nb;
			end

			IdsOp = [IdsOp idNeigh];
			%record in y
			y(counter + 1, :) = y(counter, :);
			y(counter + 1, infStatus + 1 + 3 * (1 - opNew)) = y(counter + 1, infStatus + 1 + 3 * (1 - opNew)) + 1;
			y(counter + 1, infStatus + 1 + 3 * opNew) = y(counter + 1, infStatus + 1 + 3 * opNew) - 1;

			% if more than zero people are infected
			if y(counter + 1, 2) + y(counter + 1, 5) > 0
				if people(id, 3) == 0 % the person is susceptible and their opinion switch modifies their susceptibility
					%osps(id + cfg.nPeople) = PropInfectList(People,OpNet,PhysNet,pars,id);
					osps(id + cfg.nPeople) = infectionPropensities(people, opNetRels, physNetRels, id, cfg);
				end
			end

			osps(IdsOp) = opSwitchPropensities(people, IdsOp, cfg, nInfected, ld_duration);
		else % epidemiological event
			%determine the id of individual
			id = j - cfg.nPeople;
			op = people(id, 2);
			ep = people(id, 3);

			if ep == 0 % epidemiological status is susceptible, transition to infected
				nInfected = nInfected + 1;
				people(id, 3) = ep + 1;

				%record changes in y
				y(counter + 1, :) = y(counter, :);
				y(counter + 1, (1 - op) * 3 + ep + 1) = y(counter + 1, (1 - op) * 3 + ep + 1) - 1;
				y(counter + 1, (1 - op) * 3 + ep + 2) = y(counter + 1, (1 - op) * 3 + ep + 2) + 1;

				%update the propensities for the infection of neighbours and of self
				IdsPhys = [id];
				idNeigh = FriendsPhys{id};
				idNeighOp = FriendsOp{id};

				% update number of infected and non-infected individuals for the friends of id
				nFriendsInfected = 0;
				for idN = idNeighOp
					people(idN, 6) = people(idN, 6) + 1;
					people(idN, 7) = people(idN, 7) - 1;
					if people(idN, 3) == 1
						nFriendsInfected = nFriendsInfected + 1;
					end
				end

				% log transmission number for initial growth regime
				if t(counter) < 0.01
					nTransmissions(id) = 1;
				end

				if t(counter) < 5
					iInfector = randi([1, nFriendsInfected]);
					iCounter = 1;
					for idN = idNeighOp
						if people(idN, 3) == 1
							if iCounter == iInfector
								%assert(nTransmissions(idN) > 0);
								if nTransmissions(idN) == 0
									continue;
								end
								nTransmissions(idN) = nTransmissions(idN) + 1;
								%break;
							end
							iCounter = iCounter + 1;
						end
					end
				end
				
				% update epidemiological propensity of self and people who are susceptible in the physical network retain only these that are susceptible
				idNeighSusc = idNeigh(find(people(idNeigh, 3) == 0));
				IdsPhys = [IdsPhys idNeighSusc];
				osps(IdsPhys + cfg.nPeople) = infectionPropensities(people, opNetRels, physNetRels, IdsPhys, cfg);

				%update the propensities for the opinion switch for neighbours
				osps(idNeighOp) = opSwitchPropensities(people, idNeighOp, cfg, nInfected, ld_duration);
				
			elseif ep == 1 % status is infected, transition to recovered
				nInfected = nInfected - 1;

				%update individual id
				people(id, 3) = 2;

				%update epidemiological neighbourhood
				%collect necessary ids
				IdsPhys = [id];
				idNeigh = FriendsPhys{id};
				idNeighOp = FriendsOp{id};

				% update number of infected and non-infected individuals for the
				% friends of id
				for idN = idNeighOp
					people(idN, 6) = people(idN, 6) + 1;
					people(idN, 7) = people(idN, 7) - 1;
				end

				% retain only these that are susceptible
				idNeighSusc = idNeigh(find(people(idNeigh, 3) == 0));
				IdsPhys = [IdsPhys idNeighSusc];
				osps(IdsPhys + cfg.nPeople) = infectionPropensities(people, opNetRels, physNetRels, IdsPhys, cfg);

				%record changes in y
				y(counter + 1, :) = y(counter, :);
				%recall: s - ep=0, i - ep=1, r - ep=2 but it goes 1,2,3 in y -
				%shifted by 1
				y(counter + 1, ep + 1 + 3 * (1 - op)) = y(counter + 1, ep + 1 + 3 * (1 - op)) - 1;
				y(counter + 1, ep + 2 + 3 * (1 - op)) = y(counter + 1, ep + 2 + 3 * (1 - op)) + 1;


				%update the propensities for the opinion switch for neighbours
				osps(idNeighOp) = opSwitchPropensities(people, idNeighOp, cfg, nInfected, ld_duration);
				%fprintf("I -> R\n");
			else %catchall, an error occured, since epidemiological event happened to a recovered person
				ep 
				error('SSAOpNetwork: epidemiological event has happened to a recovered individual');
			end

		end

		lds(counter + 1) = ld;

		if ~(sum(osps) > 0)
			fl = 1;
		end

		if ld == 1
			ld_duration = ld_duration + deltaT;
		end

		if (ld == 0 && nInfected / cfg.nPeople > cfg.f_s && cfg.ldSwitch)
			fprintf("lockdown initiated at prevalence %f (t = %f)\n", nInfected / cfg.nPeople, t(counter)); 
			[people, physNetRels] = initiateLockdown(people, physNetRels, cfg);
			ld = 1;
		end

		if (ld == 1 && nInfected / cfg.nPeople < cfg.f_e)
			fprintf("lockdown lifted at prevalence %f (t = %f)\n", nInfected / cfg.nPeople, t(counter));
			[people, physNetRels] = endLockdown(people, physNetRels, cfg);
			ld = 0;
			ld_duration = 0;
		end

		counter = counter + 1;
		if floor(t(counter)) == tfc
			fprintf("t: %2.2f, prev: %f\n", t(counter), nInfected / cfg.nPeople);
			tfc = tfc + 1;
		end

		if t(counter) >= tCur
			yRec = [yRec; y(counter, :)];
			tRec = [tRec; t(counter)];
			ldsRec = [ldsRec; ld];

			tCur = tCur + tRes;
		end

		%if nInfected == 0
		%	break;
		%end
	end
	fprintf("ssa end\n");
	toc

	%check whether we got to the final time, and if not 'bookend' the
	%trajectory'
	if t(counter) < cfg.T_f
		t(counter + 1) = cfg.T_f;
		y(counter + 1, :) = y(counter, :);
		yRec = [yRec; y(counter, :)];
		tRec = [tRec; cfg.T_f];
		ldsRec = [ldsRec; ldsRec(end)];
	end
end
