% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function [people, physNetRels] = initiateLockdown(people, physNetRels, cfg)
	cuts = 0;
	for i = 1:size(physNetRels, 1)
		id1 = physNetRels(i, 1);
		id2 = physNetRels(i, 2);
		cutProb = 0;
		% at least one of the two peers is of opinion A
		if people(id1, 2) + people(id2, 2) > 0
			cutProb = cfg.q;
		% both peers are of opinion B
		else 
			cutProb = (1 - cfg.alpha) * cfg.q;
		end

		rr = rand;
		if rr < cutProb
			physNetRels(i, 3) = 0;
			cuts = cuts + 1;
		end
	end

	fprintf("edges: %d, cuts: %d\n", size(physNetRels, 1), cuts);
end
