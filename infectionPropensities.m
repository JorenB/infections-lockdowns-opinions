% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function a = infectionPropensities(people, opNetRels, physNetRels, ids, cfg)
	%This function returns list a which contains propensities of an epidemiological event taking place
	% for each individual in the list ids
	%input: population table People, Opinion network OpNet, physical network
	%PhysNetMat, parameter list pars, list of individuals for whom the
	%propensities need to be recalculated

	a = zeros(1, numel(ids));
	for counter = 1:numel(ids)
		id = ids(counter);
		ep = people(id, 3);
		op = people(id, 2);

		if ep == 0 % the individual is susceptible, transition to infected
			if op == 1 % id is of opinion 'a'
				epsilon = cfg.eps_a;
			else % id is of opinion 'b'
				epsilon = cfg.eps_b;
			end

			% find infected individuals in network of id
			ids1 = [physNetRels(find(physNetRels(:, 1) == id), 2:3);  physNetRels(find(physNetRels(:, 2) == id), [1, 3])];
			ids1u = unique(sortrows(ids1), 'rows');

			infectstat = people(ids1u(:, 1), 3);
			ids1a = [ids1u'; infectstat']';
			ind = find(ids1a(:, 3) == 1);
			idsTemp = ids1a(ind, :);
			if numel(idsTemp) > 0
				a(counter) = epsilon * (sum(idsTemp(:, 2)));
			end
		elseif ep==1 % the individual is infectious, transition to recovered
			a(counter) = cfg.gamma;
		end 
	end
end
