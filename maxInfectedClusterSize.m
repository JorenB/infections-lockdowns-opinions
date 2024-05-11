function s =  maxInfectedClusterSize(people, physNetRels) 
	% augment physNetRels with infection status of each of the partners
	physNetClust = [physNetRels(:,1:2)'; people(physNetRels(:,1), 3)'; people(physNetRels(:,2), 3)']';
	%extract ids of the relationships where at least one partner is not infected
	indPhys1 = find(physNetClust(:,3) ~= 1);
	indPhys2 = find(physNetClust(:,4) ~= 1);
	%indPhys is a list of the indices of relations in physNetClust (or physNetRels) where at least one partner is not infected
	indPhys = sort(unique([indPhys1; indPhys2])); 
	% remove from physNetClust all relations where at least one partner is not infected
	physNetClust(indPhys, :) = [];

	%calculate the size of the maximum component
	% inefficient piece of code, will have to be re-written
	s = 0;
	if numel(physNetClust) > 0
		g = graph(physNetClust(:,1)',physNetClust(:,2)');
		bincell = biconncomp(g, 'OutputForm', 'cell');

		for temp_count = 1:numel(bincell)
			if numel(bincell{temp_count}) > s
				s = numel(bincell{temp_count});
			end
		end
	end
end
