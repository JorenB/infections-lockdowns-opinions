function fid = main(num)
	% this is the main file that runs an ensemble of stochastic simulations
	% using specifications provided by the configuration file number num (user
	% provided). We set out by burning in the opinion trajctory in a disease free population.
	% Once the burn in is completed, the disease is seeded and the dynamics are
	% simulated and recorded.
	% After the simulations are done running, the output is
	% collected in the form of .csv files located in the folder corresponding
	% to the scenario
	rng(0);
	

	ids = strsplit(fileread('conf/ids.txt'), ' ');
	idStr = char(ids(num));
	
	diaryFilename = ['diaries/d', idStr,'.txt'];
	yfName = ['out/y', idStr, '.csv'];
	ntName = ['out/nt', idStr, '.csv'];
	if exist(diaryFilename, 'file') ; delete(diaryFilename); end
	yfHandle = fopen(yfName, 'w');
	fclose(yfHandle);
	%ntHandle = fopen(ntName, 'w');
	%fclose(ntHandle);
	diary(diaryFilename);
	diary on;
	char(ids(num))

	cfg = cfgRead(['conf/cfg-cal', char(ids(num)), '.txt'])

	%set up the list of people
	%people_j = [id_j,op_j,infect_j,Na_j,Nb_j,Ni_j]
	people = zeros(cfg.nPeople, 7);
	people(:,1) = (1:cfg.nPeople)';

	for iTraj = 1:cfg.nTraj
		fprintf("init trajectory %d\n", iTraj);
		rng(iTraj)  % hack to ensure deterministic outcomes regardless of T_f
		% initialization with random factor present, needs to be done inside
		% the loop as to eliminate the effect of selecting specific initial
		% distribution or a specific network

		%generate network (identical networks for physical and opinion interactions)
		fprintf("SWgen\n");
		physNetAdjMat = SW(cfg.nPeople, cfg.avDegree, cfg.Beta);
		physNetRels = Adj2Rels(physNetAdjMat);
		fprintf("done\n");
		
		opNetAdjMat = physNetAdjMat;
		opNetRels = physNetRels;
		%opNetAdjMat = SW(cfg.nPeople, cfg.avDegree, cfg.Beta);
		%opNetRels = Adj2Rels(opNetAdjMat);

		% clear people data
		people(:,2:7) = zeros(cfg.nPeople, 6);

		
		people = initOpinions(people, opNetRels, cfg.initfOp_a);

		nPhysEdges = size(physNetRels, 1);
		nOpEdges = size(opNetRels, 1);

		% set physical network edge weights
		physNetRels = [physNetRels'; zeros(1, nPhysEdges) + cfg.c_phys]';
		
		% set opinion network edge weights
		opNetRels = [opNetRels'; zeros(1, nOpEdges) + cfg.c_op]';

		% burn in the opinion dynamics
		%tmpc_op = cfg.c_op;
		%cfg.c_op = 2;
		fprintf('start burn\n');
		[t, y, people] = ssaBurn(people, opNetRels, physNetRels, cfg);
		fprintf('burn completed\n');
		%cfg.c_op = tmpc_op;
		%fprintf('no burn\n');

		%seed infection
		%if cfg.initfInfected > 0
			%people = initInfect(people, cfg.initfInfected);

		people(randi([1, cfg.nPeople]), 3) = 1;
		people(randi([1, cfg.nPeople]), 3) = 1;
		people(randi([1, cfg.nPeople]), 3) = 1;
		%end

		fprintf('start disease dynamics\n');
		[t, y, lds, nTransmissions] = ssaHybrid(people, opNetRels, physNetRels, cfg);

		%clusters(end-10:end)
		%CCA(end-10:end)
		%CCB(end-10:end)
		%GCCA(end-10:end)

		yfHandle = fopen(yfName, 'a');
		fprintf(yfHandle, 't %d\n', iTraj);
		fclose(yfHandle);

		tstr = arrayfun(@(x) num2str(x, '%3.5f'), t, 'uniformoutput', false);

		yDiary = [string(tstr) y lds];
		%yTrim = yDiary(1:20:end,:);

		writematrix(yDiary, yfName, 'Delimiter', 'space', 'WriteMode', 'append');

		%writematrix(nTransmissions, ntName, 'WriteMode', 'append');

	end

	diary off;
