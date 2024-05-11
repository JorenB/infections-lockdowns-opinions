function a = opSwitchPropensities(people, ids, cfg, nInfected, ldDuration)
	%returns array a of size 1..numel(ids), where the list contains the people whose opinion propensities need to be recalculated
	% input: population table people, parameter list pars, list of ids ids,
	% total number of infected in the population infect

	p_b = cfg.p_b;
	p_a = cfg.p_a;
	c_op = cfg.c_op;
	k = cfg.k;
	theta_a = cfg.theta_a;
	theta_b = cfg.theta_b;
	nPeople = cfg.nPeople;
	C_fat = cfg.C_fat;
	C_hs = cfg.C_hs;
	q = cfg.q;
	prev_local_global = 1; % TODO in config

	a = zeros(1, numel(ids));
	for counter = 1:numel(ids)
		id = ids(counter);
		op = people(id, 2);
		Na = people(id, 4);
		Nb = people(id, 5);

		if Na + Nb > 0 % check whether the node has any neighbours            
			if op == 1% opinion a, switching to b
				% first term: switch propensity due to local opinion distribution
				% second term: switch propensity due to lockdown fatigue 
				fnb = Nb / (Na + Nb);
				a(counter) = c_op * p_b * fnb^k / (1 + theta_b * fnb^k) + C_fat * ldDuration * q;
				%a(counter) = c_op * p_b * Nb / (Na + Nb);
			else % opinion b, switching to a
				% first term: switch propensity due to local opinion distribution
				% second term: switch propensity due to health scare
				fna = Na / (Na + Nb);
				prevalence = 0;
				if prev_local_global == 0
					prevalence = people(counter, 6) / (people(counter, 6) + people(counter, 7));
				else
					prevalence = nInfected / nPeople;
				end	
				
				a(counter) = c_op * p_a * fna^k / (1 + theta_a * fna^k) + C_hs * prevalence;
				%a(counter) = c_op * p_a * Na / (Na + Nb);
			end
		else
			a(counter) = 0;
		end
	end
end
