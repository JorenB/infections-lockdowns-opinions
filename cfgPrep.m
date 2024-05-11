% this script prepares configuration files used by the main bootstrap that run
% the simulations and collect and analyze the data for the simulation where
% infection and opinion co-circulate. The following parameters are varied:
% type of network (average path length by mean of probability of edge rewiring in WS algorithm),
% social contact rate, local vs global information about prevalence, speed of reaction to prevalence

% prepare state space
clc;
clear variables;
close all;
format long;

folderStr = 'conf';
mkdir(folderStr)
fileStr = 'cfg-';

% lockdown switch
%ldArr = [0, 1];
ldArr = [0];

% stringency rates
%qArr = [0.2, 0.5, 0.9];
%qArr = [0.1, 0.3, 0.5, 0.7, 0.9];
%qArr = [0, 0.25, 0.5, 0.75, 1];
qArr = [0.375, 0.625];

% lockdown resistance
%alphaArr = [0.15, 0.3, 0.45, 0.6, 0.75];
%alphaArr = [0, 0.5, 1];
alphaArr = [0];

% q_alpha_Arr = [[0.0, 0.0]; [0.1, 0.0]; [0.3, 0.0]; [0.5, 0.0]; [0.7, 0.0]; [0.9, 0.0]; [1.0, 0.0]; [0.1, 0.9]; [0.3, 0.9]; [0.5, 0.9]; [0.7, 0.9]; [0.9, 0.9]; [1.0, 0.15]; [1.0, 0.3]; [1.0, 0.45]; [1.0, 0.6]; [1.0, 0.75]; [1.0, 0.9]; [1.0, 1.0]]';

%deg_cr_eps_ab_chs_clfArr = [[14, 7.0, 0.02, 0.035, 2, 0.05]; [38, 2.5, 0.015, 0.03, 2, 0.05]]';
deg_cr_eps_ab_chs_clfArr = [[14, 7.0, 0.02, 0.035, 2, 0.05]]';

fs_Arr = [0.005, 0.01, 0.02, 0.035, 0.05];
%fs_Arr = [0.0425];

ids = [];
lookup = [];

nPeople = 10000;

c = 1;
for deg_cr_eps_ab_chs_clf = deg_cr_eps_ab_chs_clfArr
	for q = qArr
		for alpha = alphaArr
	%for qa = q_alpha_Arr
		for fs = fs_Arr
			cfg = jsondecode(fileread('cfgTemplate.txt'));
			cfg.nPeople = nPeople;
			cfg.ldSwitch = 1;
			cfg.avDegree = deg_cr_eps_ab_chs_clf(1);
			cfg.c_phys = deg_cr_eps_ab_chs_clf(2);
			cfg.eps_a = deg_cr_eps_ab_chs_clf(3);
			cfg.eps_b = deg_cr_eps_ab_chs_clf(4);
			cfg.C_hs = deg_cr_eps_ab_chs_clf(5);
			cfg.C_fat = deg_cr_eps_ab_chs_clf(6);
			cfg.initfOp_a = 0.5;
			cfg.tBurn = 100;
			cfg.T_f = 50;
			% cfg.q = qa(1);
			% cfg.alpha = qa(2);
			cfg.q = q;
			cfg.alpha = alpha;
			cfg.f_s = fs

			id = ['-deg', num2str(cfg.avDegree), '-N', num2str(nPeople), '-q', num2str(cfg.q), '-alpha', num2str(cfg.alpha), '-fs', num2str(cfg.f_s)];
			ids = [ids, ' ', id];
			lookup = [lookup, id, ' ', num2str(c), '\n'];
			c = c + 1;

			fname = ['cfg-cal', id, '.txt'];

			fid = fopen(['conf/', fname], 'w');
			fprintf(fid, jsonencode(cfg));
			fclose(fid);
		end
	end
	end
end
%for ld = ldArr
%	for q = qArr
%		for alpha = alphaArr
%			for eps_a = eps_aArr
%				cfg = jsondecode(fileread('cfgTemplate.txt'));
%				cfg.nPeople = nPeople;
%				cfg.ldSwitch = ld;
%				cfg.q = q;
%				cfg.alpha = alpha;
%				cfg.eps_a = eps_a;
%				cfg.eps_b = 2 * eps_a;
%
%				id = ['-ld', num2str(ld), '-q', num2str(q), '-alpha', num2str(alpha), '-N', num2str(nPeople), '-ea', num2str(eps_a), '-eb', num2str(cfg.eps_b)];
%				ids = [ids, ' ', id];
%				lookup = [lookup, id, ' ', num2str(c), '\n'];
%				c = c + 1;
%
%				fname = ['cfg', id, '.txt'];
%
%				fid = fopen(['conf/', fname], 'w');
%				fprintf(fid, jsonencode(cfg));
%				fclose(fid);
%			end
%		end
%	end
%end


disp(lookup);
fid = fopen('conf/ids.txt', 'w');
fprintf(fid, ids(2:end));
fclose(fid);

fid = fopen('conf/lookup.txt', 'w');
fprintf(fid, lookup);
fclose(fid);
