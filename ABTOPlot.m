function res=ABTOPlot(pars)
% this function outputs the time series for the coupled model of opinion 
% competition with infection spread (SIR) without adaptive behavior or
% assortativity. The data used is recorded in the in Scenarios
% folder and will be presented in the shape of arrays and figures
close all;
res=0;
str='Scenarios/Specs_';

% so far no containers are initiated
% format strings definition
copArr=[0.2;1;2];

PrevFlarr=[0,1];

marr=[1,25,75];

betaLow=-4;
betaUp=-1;
betaNum=6;
betaArr=10.^linspace(betaLow,betaUp,betaNum);
betaArr=[0,betaArr];
% cop x NaFl']=[0.2,0; 1,0; 2,0; 0.2,1; 1,1; 2,1]
%define colors
pink=[239,38,72]/255;
Lpink=[255,158,168]/255;
teal=[90,138,145]/255;
Lteal=[169,255,228]/255;
mustard=[251,200,0]/255;
Lmustard=[252,255,172]/255;

colors=[pink;teal;mustard];
Lcolors=[Lpink;Lteal;Lmustard];

counter=1;

% initialize containers for the long term dynamics
% opinion dynamics
% density of Na
NaPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
NaPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));

% the size of the largest biconnected B-component
BCPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
BCPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));

%local clustering coefficient of Na network
LCPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
LCPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));

%global clustering coefficient of Na network
GCPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
GCPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));

%infection dynamics

%local clustering coefficient of the infected network
LCCInfectPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
LCCInfectPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));
%final size
FSPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
FSPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));
% the size of the largest biconnected component in the infected network
BCInfectPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
BCInfectPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));
%prevalence of infectious individuals
PrevPrevFl0=zeros(numel(copArr),betaNum+1,numel(marr));
PrevPrevFl1=zeros(numel(copArr),betaNum+1,numel(marr));

for c0=1:1:3 % counter of m
    for c1=1:1:3 %contact rate counter
        for c2=1:2 %PrevFl counter
            for c3=1:(betaNum+1) %Network counter
                %compose file names
                num=c3+((c0-1)*4+(c1-1)*2+c2-1)*(betaNum+1);
                tstr=[str,num2str(num),'/tdiscr.csv'];
                %opinion strings
                Nastr=[str,num2str(num),'/Na.csv'];
                CCAstr=[str,num2str(num),'/CCA.csv'];
                GCCAstr=[str,num2str(num),'/GCCA.csv'];
                Clusterstr=[str,num2str(num),'/Cluster.csv'];
                %infection strings
                LCCInfectstr=[str,num2str(num),'/CCInfect.csv'];
                FSstr=[str,num2str(num),'/FS.csv'];
                BCInfectstr=[str,num2str(num),'/InfectCluster.csv'];
                PrevPrevstr=[str,num2str(num),'/Infect.csv'];
                %network strings
                GCCstr=[str,num2str(num),'/GCCNetwork.csv'];
                ShortestPstr=[str,num2str(num),'/ShortestP.csv'];

                %load files
                t=readmatrix(tstr);
                % opinion dynamics
                NA=readmatrix(Nastr);
                CCA=readmatrix(CCAstr);
                GCCA=readmatrix(GCCAstr);
                Cluster=readmatrix(Clusterstr);
                % infection dynamics
                LCCInfect=readmatrix(LCCInfectstr);
                FSPrev=readmatrix(FSstr);
                BCInfect=readmatrix(BCInfectstr);
                PrevPrev=readmatrix(PrevPrevstr);
                [num_traj,num_points]=size(PrevPrev);
                Prob_ext=zeros(num_traj,num_points);
                Peak=zeros(num_traj,num_points);


                GCC=readmatrix(GCCstr);
                ShortestP=readmatrix(ShortestPstr);
    
                % calculate medians and 10-th and 90-th quantile
                % opinion dynamics
                NA_10q=quantile(NA,0.1);
                NA_50q=quantile(NA,0.5);
                NA_90q=quantile(NA,0.9);
                
                CCA_10q=quantile(CCA,0.05);
                CCA_50q=quantile(CCA,0.5);
                CCA_90q=quantile(CCA,0.95);
    
                GCCA_10q=quantile(GCCA,0.05);
                GCCA_50q=quantile(GCCA,0.5);
                GCCA_90q=quantile(GCCA,0.95);
    
                Cluster_10q=quantile(Cluster,0.1);
                Cluster_50q=quantile(Cluster,0.5);
                Cluster_90q=quantile(Cluster,0.9);

                %infection dynamics
                LCCInfect_10q=quantile(LCCInfect,0.1);
                LCCInfect_50q=quantile(LCCInfect,0.5);
                LCCInfect_90q=quantile(LCCInfect,0.9);

                FSPrev_10q=quantile(FSPrev,0.1);
                FSPrev_50q=quantile(FSPrev,0.5);
                FSPrev_90q=quantile(FSPrev,0.9);

                BCInfect_10q=quantile(BCInfect,0.1);
                BCInfect_50q=quantile(BCInfect,0.5);
                BCInfect_90q=quantile(BCInfect,0.9);

                PrevPrev_10q=quantile(PrevPrev,0.1);
                PrevPrev_50q=quantile(PrevPrev,0.5);
                PrevPrev_90q=quantile(PrevPrev,0.9);

                    
%                 % plot depending on c1, c2, c3
%                 %density of N_{a}                
%                 figure((c3+c0*(betaNum+1))*9+1);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[100*NA_10q, fliplr(100*NA_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,100*NA_10q,'Color',colors(c1,:),'LineWidth',1)
%       
%                 plot(t,100*NA_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==2
%                     hNa(c3,c1)=plot(t,100*NA_50q,'Color',colors(c1,:),'LineWidth',4);
%                 else
%                     plot(t,100*NA_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('Density of $$N_{a}$$, $$n_{a}$$ ($$\%$$)','Interpreter','latex');
%                 end
%                 if c2==1
%                     title('Random placement','Interpreter','latex');
%                 else
%                     title('Largest degree placement','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
%                 
%                 %local cluster coefficient of N_{a} network
%                 figure((c3+c0*(betaNum+1))*9+2);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[CCA_10q, fliplr(CCA_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,CCA_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hCCA(c3,c1)=plot(t,CCA_50q,'Color',colors(c1,:),'LineWidth',4);
%                 else
%                     plot(t,CCA_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('LCC of $$N_{a}$$','Interpreter','latex');
%                 end
%                 plot(t,CCA_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 
%                 if c2==1
%                     title('Random placement','Interpreter','latex');
%                 else
%                     title('Largest degree placement','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
%     
%                 %global cluster coefficient of N_{a} network
%                 figure((c3+c0*(betaNum+1))*9+3);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[GCCA_10q, fliplr(GCCA_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,GCCA_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hGCCA(c3,c1)=plot(t,GCCA_50q,'Color',colors(c1,:),'LineWidth',4);
%                 else
%                     plot(t,GCCA_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('GCC of $$N_{a}$$','Interpreter','latex');
%                 end
%                 plot(t,GCCA_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 
%                 if c2==1
%                     title('Random placement','Interpreter','latex');
%                 else
%                     title('Largest degree placement','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
%     
%     
%                 %the size of the largest biconnected component of N_{b}
%                 %network
%                 figure((c3+c0*(betaNum+1))*9+4);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[Cluster_10q, fliplr(Cluster_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,Cluster_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hCluster(c3,c1)=plot(t,Cluster_50q,'Color',colors(c1,:),'LineWidth',4);
%                 else
%                     plot(t,Cluster_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel({'Size of the largest';'biconnected component of $$N_{b}$$'},'Interpreter','latex');
%                 end
%                 plot(t,Cluster_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Random placement','Interpreter','latex');
%                 else
%                     title('Largest degree placement','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
%             
%                 %% infection dynamics
%                 % local clustering coefficient
%                 figure((c3+c0*(betaNum+1))*9+5);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[LCCInfect_10q, fliplr(LCCInfect_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,LCCInfect_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hLCCInfect(c3+c0*(betaNum+1),c1)=plot(t,LCCInfect_50q,'Color',colors(c1,:),'LineWidth',4); 
%                 else
%                     plot(t,LCCInfect_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel({'LCC of';'infected network'},'Interpreter','latex');
%                 end
%                 plot(t,LCCInfect_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Local information','Interpreter','latex');
%                 else
%                     title('Global information','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
%                 
%                 % final size
%                 figure((c3+c0*(betaNum+1))*9+6);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[FSPrev_10q, fliplr(FSPrev_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,FSPrev_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hFSPrev(c3+c0*(betaNum+1),c1)=plot(t,FSPrev_50q,'Color',colors(c1,:),'LineWidth',4); 
%                 else
%                     plot(t,FSPrev_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('Final size','Interpreter','latex');
%                 end
%                 plot(t,FSPrev_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Local information','Interpreter','latex');
%                 else
%                     title('Global information','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
% 
%                 % the size of the largest biconnected component in
%                 % infectious network
%                 figure((c3+c0*(betaNum+1))*9+7);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[BCInfect_10q, fliplr(BCInfect_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,BCInfect_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hBCInfect(c3+c0*(betaNum+1),c1)=plot(t,BCInfect_50q,'Color',colors(c1,:),'LineWidth',4); 
%                 else
%                     plot(t,BCInfect_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel({'Size of the largest';'biconnected component';' of infected network'},'Interpreter','latex');
%                 end
%                 plot(t,BCInfect_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Local information','Interpreter','latex');
%                 else
%                     title('Global information','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
% 
%                 
%                 figure((c3+c0*(betaNum+1))*9+8);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[PrevPrev_10q, fliplr(PrevPrev_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,PrevPrev_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hPrevPrev(c3+c0*(betaNum+1),c1)=plot(t,PrevPrev_50q,'Color',colors(c1,:),'LineWidth',4); 
%                 else
%                     plot(t,PrevPrev_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('Prevalence of infectious individuals','Interpreter','latex');
%                 end
%                 plot(t,PrevPrev_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Local information','Interpreter','latex');
%                 else
%                     title('Global information','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
% 
                %% calculate probability of extinction for each trajectory
                [num_traj,num_points]=size(PrevPrev);

                Prob_ext=zeros(num_traj,num_points);
                PeakSize=zeros(1,num_traj);
                PeakTime=zeros(1,num_traj);
                for tcounter0=1:1:num_traj
                    ind=find(PrevPrev(tcounter0,:)==0,1);
                    onet=ones(1,num_points-ind+1);
                    onet=padarray(onet,[0 ind-1],0,'pre');
                    Prob_ext(tcounter0,:)=onet;
                    clear onet;
                    % calculate peak and time of its occurence
                    ind=find(PrevPrev(tcounter0,:)==max(PrevPrev(tcounter0,:)),1);
                    PeakSize(tcounter0)=PrevPrev(tcounter0,ind);
                    PeakTime(tcounter0)=t(ind);
                end

                % calculate median and quantile for probability of
                % extinction
                Prob_ext_10q=quantile(Prob_ext,0.1);
                Prob_ext_50q=quantile(Prob_ext,0.5);
                Prob_ext_90q=quantile(Prob_ext,0.9);
%                 
%                 figure((c3+c0*(betaNum+1))*9+9);
%                 subplot(1,2,c2)
%                 h=fill([t,fliplr(t)],[Prob_ext_10q, fliplr(Prob_ext_90q)],Lcolors(c1,:));
%                 set(h,'facealpha',.5)
%                 hold on;
%                 plot(t,Prob_ext_10q,'Color',colors(c1,:),'LineWidth',1)
%                 if c2==2
%                     hProb_ext(c3+c0*(betaNum+1),c1)=plot(t,Prob_ext_50q,'Color',colors(c1,:),'LineWidth',4); 
%                 else
%                     plot(t,Prob_ext_50q,'Color',colors(c1,:),'LineWidth',4);
%                     ylabel('Probability of infection extinctions','Interpreter','latex');
%                 end
%                 plot(t,Prob_ext_90q,'Color',colors(c1,:),'LineWidth',1)
%                 xlabel('Time, weeks','Interpreter','latex');
%                 if c2==1
%                     title('Local information','Interpreter','latex');
%                 else
%                     title('Global information','Interpreter','latex');
%                 end
%                 set(gca,'FontSize',25);
                if (c3==2) || (c3==6) % output time series
                    % opinion dynamics
                    % density of N_{a}                
                    figure((c3+c0*(betaNum+1))*9+1);
                    
                    subplot(1,2,c2)
                    h=fill([t,fliplr(t)],[NA_10q, fliplr(NA_90q)],Lcolors(c1,:));
                    set(h,'facealpha',.5);
                    hold on;
                    if c2==1 && c1==1 % write annotation box
                        strann=['$$\beta=',num2str(betaArr(c3),2),'$$ and $$m=',num2str(marr(c0)),'$$'];
                        annotation('textbox',[0.3,0,1,1],'String',strann,'FitBoxToText','on','FontSize',25,'interpreter','latex','FontWeight','bold','EdgeColor','none');
                    end    
                    
                    plot(t,NA_10q,'Color',colors(c1,:),'LineWidth',1);
                    plot(t,NA_90q,'Color',colors(c1,:),'LineWidth',1)
          
                    xlabel('Time, weeks','Interpreter','latex');
                    if c2==2
                        hNa(c3+c0*(betaNum+1),c1)=plot(t,NA_50q,'Color',colors(c1,:),'LineWidth',4);
                    else
                        plot(t,NA_50q,'Color',colors(c1,:),'LineWidth',4);
                        ylabel('Density of $$N_{a}$$, $$n_{a}$$','Interpreter','latex');
                    end
                    if c2==1
                        title('Local information','Interpreter','latex');
                    else
                        title('Global information','Interpreter','latex');
                    end
                    set(gca,'FontSize',25);

                    % local cluster coefficient of N_{a} network
                    figure((c3+c0*(betaNum+1))*9+2);
                    subplot(1,2,c2)
                    h=fill([t,fliplr(t)],[CCA_10q, fliplr(CCA_90q)],Lcolors(c1,:));
                    set(h,'facealpha',.5)
                    hold on;
                    plot(t,CCA_10q,'Color',colors(c1,:),'LineWidth',1)
                    if c2==2
                        hCCA(c3+c0*(betaNum+1),c1)=plot(t,CCA_50q,'Color',colors(c1,:),'LineWidth',4);
                    else
                        plot(t,CCA_50q,'Color',colors(c1,:),'LineWidth',4);
                        ylabel('LCC of $$N_{a}$$','Interpreter','latex');
                    end

                    if c2==1 && c1==1 % write annotation box
                        strann=['$$\beta=',num2str(betaArr(c3),2),'$$ and $$m=',num2str(marr(c0)),'$$'];
                        annotation('textbox',[0.3,0,1,1],'String',strann,'FitBoxToText','on','FontSize',25,'interpreter','latex','FontWeight','bold','EdgeColor','none');
                    end
                    plot(t,CCA_90q,'Color',colors(c1,:),'LineWidth',1)
                    xlabel('Time, weeks','Interpreter','latex');
                    
                    if c2==1
                        title('Local information','Interpreter','latex');
                    else
                        title('Global information','Interpreter','latex');
                    end
                    set(gca,'FontSize',25);

                    % infection dynamics
                    % prevalence
                    figure((c3+c0*(betaNum+1))*9+8);
                    subplot(1,2,c2)
                    h=fill([t,fliplr(t)],[PrevPrev_10q, fliplr(PrevPrev_90q)],Lcolors(c1,:));
                    set(h,'facealpha',.5)
                    hold on;
                    plot(t,PrevPrev_10q,'Color',colors(c1,:),'LineWidth',1)
                    if c2==2
                        hPrevPrev(c3+c0*(betaNum+1),c1)=plot(t,PrevPrev_50q,'Color',colors(c1,:),'LineWidth',4); 
                    else
                        plot(t,PrevPrev_50q,'Color',colors(c1,:),'LineWidth',4);
                        ylabel('Prevalence of infectious individuals','Interpreter','latex');
                    end
                    if c2==1 && c1==1 % write annotation box
                        strann=['$$\beta=',num2str(betaArr(c3),2),'$$ and $$m=',num2str(marr(c0)),'$$'];
                        annotation('textbox',[0.3,0,1,1],'String',strann,'FitBoxToText','on','FontSize',25,'interpreter','latex','FontWeight','bold','EdgeColor','none');
                    end
                    plot(t,PrevPrev_90q,'Color',colors(c1,:),'LineWidth',1)
                    xlabel('Time, weeks','Interpreter','latex');
                    if c2==1
                        title('Local information','Interpreter','latex');
                    else
                        title('Global information','Interpreter','latex');
                    end
                    set(gca,'FontSize',25);

                    % probability of extinction
                    figure((c3+c0*(betaNum+1))*9+9);
                    subplot(1,2,c2)
                    h=fill([t,fliplr(t)],[Prob_ext_10q, fliplr(Prob_ext_90q)],Lcolors(c1,:));
                    set(h,'facealpha',.5)
                    hold on;
                    plot(t,Prob_ext_10q,'Color',colors(c1,:),'LineWidth',1)
                    if c2==2
                        hProb_ext(c3+c0*(betaNum+1),c1)=plot(t,Prob_ext_50q,'Color',colors(c1,:),'LineWidth',4); 
                    else
                        plot(t,Prob_ext_50q,'Color',colors(c1,:),'LineWidth',4);
                        ylabel('Probability of infection extinctions','Interpreter','latex');
                    end
                    if c2==1 && c1==1 % write annotation box
                        strann=['$$\beta=',num2str(betaArr(c3),2),'$$ and $$m=',num2str(marr(c0)),'$$'];
                        annotation('textbox',[0.3,0,1,1],'String',strann,'FitBoxToText','on','FontSize',25,'interpreter','latex','FontWeight','bold','EdgeColor','none');
                    end
                    plot(t,Prob_ext_90q,'Color',colors(c1,:),'LineWidth',1)
                    xlabel('Time, weeks','Interpreter','latex');
                    if c2==1
                        title('Local information','Interpreter','latex');
                    else
                        title('Global information','Interpreter','latex');
                    end
                    set(gca,'FontSize',25);
                end    
                
                % collect end of the run statistics for the infection
                if c2==1
                    % opinion 4 measurements
                    NaNaFl0(c1,c3+(c0-1)*(betaNum+1))=NA_50q(end); 
                    BCNaFl0(c1,c3+(c0-1)*(betaNum+1))=Cluster_50q(end);
                    LCNaFl0(c1,c3+(c0-1)*(betaNum+1))=CCA_50q(end);
                    GCNaFl0(c1,c3+(c0-1)*(betaNum+1))=GCCA_50q(end);
                    %infection 7 measurements
                    LCCInfectFl0(c1,c3+(c0-1)*(betaNum+1))=LCCInfect_50q(end); 
                    FSPrevFl0(c1,c3+(c0-1)*(betaNum+1))=FSPrev_50q(end);
                    BCInfect0(c1,c3+(c0-1)*(betaNum+1))=BCInfect_50q(end);
                    PrevPrev0(c1,c3+(c0-1)*(betaNum+1))=PrevPrev_50q(end);
                    Prob_ext0(c1,c3+(c0-1)*(betaNum+1))=Prob_ext_50q(end);
                    PeakSize0(c1,c3+(c0-1)*(betaNum+1))=quantile(PeakSize,0.5);
                    PeakTime0(c1,c3+(c0-1)*(betaNum+1))=quantile(PeakTime,0.5);
                else
                    %opinion 4 measurements
                    NaNaFl1(c1,c3+(c0-1)*(betaNum+1))=NA_50q(end);
                    BCNaFl1(c1,c3+(c0-1)*(betaNum+1))=Cluster_50q(end);
                    LCNaFl1(c1,c3+(c0-1)*(betaNum+1))=CCA_50q(end);
                    GCNaFl1(c1,c3+(c0-1)*(betaNum+1))=GCCA_50q(end);
                    % infection 7 measurements
                    LCCInfectFl1(c1,c3+(c0-1)*(betaNum+1))=LCCInfect_50q(end);
                    FSPrevFl1(c1,c3+(c0-1)*(betaNum+1))=FSPrev_50q(end);
                    BCInfect1(c1,c3+(c0-1)*(betaNum+1))=BCInfect_50q(end);
                    PrevPrev1(c1,c3+(c0-1)*(betaNum+1))=PrevPrev_50q(end);
                    Prob_ext1(c1,c3+(c0-1)*(betaNum+1))=Prob_ext_50q(end);
                    PeakSize1(c1,c3+(c0-1)*(betaNum+1))=quantile(PeakSize,0.5);
                    PeakTime1(c1,c3+(c0-1)*(betaNum+1))=quantile(PeakTime,0.5);
                end
                                
                counter=counter+1;
            end
        end
    end
end

% %output legends
for c0=1:1:3
    for c3=1:(betaNum+1)
        if c3==2 || c3==6
            if c0==1
                 figure((c3+c0*(betaNum+1))*9+1);
                 legend(hNa(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+2);
                 legend(hCCA(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+8);
                 legend(hPrevPrev(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+9);
                 legend(hProb_ext(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
            else % temporary
                 figure((c3+c0*(betaNum+1))*9+1);
                 legend(hNa(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+2);
                 legend(hCCA(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+8);
                 legend(hPrevPrev(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
                 figure((c3+c0*(betaNum+1))*9+9);
                 legend(hProb_ext(c3+c0*(betaNum+1),:),'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
            end
        end
    end
end

%%
%output summary of the median outputs at the end of simulations

% 4+ 7 figures
for c1=1:1:numel(copArr)
    for c0=1:1:3 % rows
        figc1=401;
        figure(figc1);
        subplot(3,2,(c0-1)*2+1);
        if c0==1
            hsumNa(c1,:)=plot(log(betaArr)/log(10),NaNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        else
            plot(log(betaArr)/log(10),NaNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        end
        ylim([0.5,0.56]);
        
        ylabel({'Density of $$N_{a}$$';'individuals, $$n_{a}$$';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
        if c0==1
            title('Local information','Interpreter','latex');
        end

        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        set(gca,'FontSize',25)
        subplot(3,2,(c0-1)*2+2);
        plot(log(betaArr)/log(10),NaNaFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        ylim([0.5,0.56]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        if c0==1
            title('Global information','Interpreter','latex');
        end
        set(gca,'FontSize',25);
    
%         figc1=402;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         if c0==1
%             hsumBCon(c1,:)=plot(log(betaArr)/log(10),BCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         else
%             plot(log(betaArr)/log(10),BCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         end
%         ylim([0,500]);
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'Size of the ';'largest biconnected';'component of $$N_{b}$$';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25)
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),BCNaFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         ylim([0,500]);
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

        figc1=403;
        figure(figc1);
        subplot(3,2,(c0-1)*2+1);
        if c0==1
            hsumLccNa(c1,:)=plot(log(betaArr)/log(10),LCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        else
            plot(log(betaArr)/log(10),LCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        end
        ylim([0.25,0.38]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        ylabel({'LCC of $$N_{a}$$';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
        if c0==1
            title('Local information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

        subplot(3,2,(c0-1)*2+2);
        plot(log(betaArr)/log(10),LCNaFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        ylim([0.25,0.38]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        if c0==1
            title('Global information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

%         figc1=404;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         if c0==1
%             hsumGccNa(c1,:)=plot(log(betaArr)/log(10),GCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         else
%             plot(log(betaArr)/log(10),GCNaFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         end
%         ylim([0.5,0.71]);
%         xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         ylabel({'GCC of $$N_{a}$$';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25)
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),GCNaFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         ylim([0.5,0.71]);
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

%         figc1=405;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         if c0==1
%             hsumLccInfect(c1,:)=plot(log(betaArr)/log(10),LCCInfectFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         else
%             plot(log(betaArr)/log(10),LCCInfectFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         end
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'LCC of';'infected network';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);
% 
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),LCCInfectFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

        figc1=406;
        figure(figc1);
        subplot(3,2,(c0-1)*2+1);
        if c0==1
            hsumFS(c1,:)=semilogy(log(betaArr)/log(10),FSPrevFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        else
            semilogy(log(betaArr)/log(10),FSPrevFl0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        end
        ylim([0,1]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        ylabel({'$$\log$$ Final size';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
        if c0==1
            title('Local information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

        subplot(3,2,(c0-1)*2+2);
        plot(log(betaArr)/log(10),FSPrevFl1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        ylim([0.01,0.02]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        if c0==1
            title('Global information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

%         figc1=407;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         if c0==1
%             hsumBCInfect(c1,:)=plot(log(betaArr)/log(10),BCInfect0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         else
%             plot(log(betaArr)/log(10),BCInfect0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         end
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'Size of the largest';'biconnected component';' of infected network';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);
% 
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),BCInfect1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

%         figc1=408;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         if c0==1
%             hsumPrev(c1,:)=plot(log(betaArr)/log(10),PrevPrev0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         else
%             plot(log(betaArr)/log(10),PrevPrev0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         end
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'Prevalence of infectious individuals';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),PrevPrev1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

%         figc1=409;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         hsumProbExt(c1,:)=plot(log(betaArr)/log(10),Prob_ext0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'Probability of infection extinction';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);
% 
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),Prob_ext1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);

        figc1=410;
        figure(figc1);
        subplot(3,2,(c0-1)*2+1);
        hsumPS(c1,:)=semilogy(log(betaArr)/log(10),PeakSize0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        ylim([0,0.079]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        ylabel({'$$\log$$ Peak size of infection';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
        if c0==1
            title('Local information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

        subplot(3,2,(c0-1)*2+2);
        semilogy(log(betaArr)/log(10),PeakSize1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
        ylim([0,0.079]);
        if c0==3
            xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
        end 
        if c0==1
            title('Global information','Interpreter','latex');
        end
        set(gca,'FontSize',25);

%         figc1=411;
%         figure(figc1);
%         subplot(3,2,(c0-1)*2+1);
%         hsumPT(c1,:)=plot(log(betaArr)/log(10),PeakTime0(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         ylabel({'Peak time of infection';' ';'Sensitivity of reaction'; ['$$m=',num2str(num2str(marr(c0))),'$$']},'Interpreter','latex');
%         if c0==1
%             title('Local information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);
% 
%         subplot(3,2,(c0-1)*2+2);
%         plot(log(betaArr)/log(10),PeakTime1(c1,((betaNum+1)+(c0-1)*(betaNum+1)-7+1):1:((betaNum+1)+(c0-1)*(betaNum+1))),'color',colors(c1,:),'LineWidth',4);hold on;
%         if c0==3
%             xlabel('$$\log_{10}\ \beta$$','interpreter','latex');
%         end 
%         if c0==1
%             title('Global information','Interpreter','latex');
%         end
%         set(gca,'FontSize',25);
     end
end

% set up the legends for the summary figures

figure(401);
legend(hsumNa,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(402);
% legend(hsumBCon,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
figure(403);
legend(hsumLccNa,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(404);
% legend(hsumGccNa,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(405);
% legend(hsumLccInfect,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
figure(406);
legend(hsumFS,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(407);
% legend(hsumBCInfect,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(408);
% legend(hsumPrev,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(409);
% legend(hsumProbExt,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
figure(410);
legend(hsumPS,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
% figure(411);
% legend(hsumPT,'$$c_{op}=0.2$$','$$c_{op}=1$$','$$c_{op}=2$$','interpreter','latex','Location','southoutside','Orientation','horizontal');
end
