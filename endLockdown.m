function [people, physNetRels] = endLockdown(people, physNetRels, cfg)
	physNetRels(:, 3) = cfg.c_phys;
end
