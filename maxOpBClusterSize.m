% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function s =  maxOpBClusterSize(people, opNetRels) 
	% augment opNetRels with opinions of each of the partners
	opNetClust = [opNetRels(:,1:2)'; people(opNetRels(:,1),2)'; people(opNetRels(:,2),2)']';
	%extract ids of the relationships where at least one partner believes 'a'
	indOps1 = find(opNetClust(:,3) == 1);
	indOps2 = find(opNetClust(:,4) == 1);
	%indOps is a list of the indices of relations in opNetClust (or opNetRels) where at least one partner believes 'a'
	indOps = sort(unique([indOps1; indOps2])); 
	% remove from opNetClust all relations where at least one partner believes 'a'
	opNetClust(indOps, :) = [];

	%calculate the size of the maximum component
	% inefficient piece of code, will have to be re-written
	s = 0;
	if numel(opNetClust) > 0
		g = graph(opNetClust(:,1)',opNetClust(:,2)');
		bincell = biconncomp(g, 'OutputForm', 'cell');

		for temp_count = 1:numel(bincell)
			if numel(bincell{temp_count}) > s
				s = numel(bincell{temp_count});
			end
		end
	end
end
