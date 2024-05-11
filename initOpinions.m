function people = initOpinions(people, opNetRels, initfOp_a)
	nPeople = size(people, 1);
	nOp_a = floor(nPeople * initfOp_a);

	ids = datasample(1:nPeople, nOp_a, 'Replace', false);
	people(ids, 2) = 1;
