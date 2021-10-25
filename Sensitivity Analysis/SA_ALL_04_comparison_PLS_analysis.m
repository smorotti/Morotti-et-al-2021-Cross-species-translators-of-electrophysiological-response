% This file performs the linear regression analysis and compares results
% in mouse, rabbit, and human.

close all
clear
clc

%% Output selection
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp

output_selection = [5 7 11 14]; % 4 outputs
% APD90	APD50	CaTamp	CaTtau

N_outputs = length(output_selection);

if N_outputs==1 || N_outputs==2
    sp1 = 2; sp2 = 1;
    max_panels = 2;
end
if N_outputs==3 || N_outputs==4
    sp1 = 2; sp2 = 2;
    max_panels = 4;
end
if N_outputs>4
    sp1 = 2; sp2 = 3;
    max_panels = 6;
end

N_figures = ceil(N_outputs/max_panels);

%% Loading
load parameter_matrix_1000

disp('Control, 1-Hz pacing')

load outputs_matrix_1000_0p1_300s_mouse
all_outputs_1 = all_outputs;

load outputs_matrix_1000_0p1_300s_rabbit
all_outputs_2 = all_outputs;

load outputs_matrix_1000_0p1_300s_human
all_outputs_3 = all_outputs;

disp('*****************')

color_1 = [0.470588235294118 0.670588235294118 0.188235294117647]; % GREEN
color_2 = [0.850980392156863 0.329411764705882 0.101960784313725]; % ORANGE
color_3 = [0 0.450980392156863 0.741176470588235]; % BLUE

%% Load parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
% and array 'parameter_names'
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak

par_selection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]; % all
%par_selection = [1 2 3 4 5 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23]; % MOUSE
%par_selection = [1 2 3 4 5 6 7 8 12 13 14 15 16 17 18 19 20 21 22 23]; % RABBIT
%par_selection = [1 2 3 4 5 6 7 8 12 13 15 16 17 18 19 20 21 22 23]; % HUMAN
    %par_selection = [1 2 3 4 5 8 12 13 15 16 17 18 19 20 21 22 23]; % COMMON
    
all_parameters = all_parameters(:,par_selection);
parameter_names = parameter_names(par_selection);

[N_trials N_pars] = size(all_parameters);

%% Load outputs
% load matrix all_outputs (columns: N outputs, rows: N trials)
% and array 'output_names' (and 'output_units')

% Select outputs
%all_outputs = all_outputs(:,output_selection);
output_names = output_names(output_selection);
output_units = output_units(output_selection);

N_outputs = length(output_names);

%% MOUSE 

all_outputs = all_outputs_1(:,output_selection);

% Outputs - population level
%all_outputs_mean = mean(all_outputs);
%all_outputs_std_dev = std(all_outputs);

% Check basic properties, and define X and Y
% Exclude cells with outputs = 0
%good_trials = all_outputs(:,1)>0;
    good_trials = (all_outputs(:,1)>0) .* (all_outputs(:,2)<400); % APD<400ms
good_trials_log = logical(good_trials);
% Good count: number of good trials
good_count = sum(good_trials_log);
% Good_parameters: array with parameters from good trials only
good_parameters = all_parameters(good_trials_log,:);
X = log(good_parameters);
% Good_outputs: array with parameters from good trials only
good_outputs = all_outputs(good_trials_log,:);
Y = log(good_outputs);

%all_outputs_mean = mean(good_outputs);
%all_outputs_std_dev = std(good_outputs);

% Istograms
dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_1)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% % PLS - svds algorithm (2010)
% [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,Yjack,R2x,R2y,...
%           RESSy,PRESSy,Q2,r2y_random,rv_random,Yhat4Press,Yhat4Ress] = ...
%           PLS_jack_svds(X,Y,rank(X));

% N_pars, N_outputs, N_trials, good_count

% Calculate agreement of values predicted by regression (Yhat = Bpls*X) 
% with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

avg_R2 = mean(R2each);

% R2
figure(100),set(gcf,'color','w','Position',[50,100,1500,750])
subplot(3,1,1)
bar(R2each,'FaceColor',color_1)
set(gca,'box','off','tickdir','out','fontsize',12)
set(gca,'XTick',1:N_outputs)
set(gca,'XTickLabel',output_names)
set(gca,'XLim',[0 N_outputs+1])
ylim([0 1])
rotateXLabels( gca(), 90)
%title('R^2 values')
title(['R^2 values (mean = ',num2str(avg_R2,3),')'])

% dex1 = 1;
% for figdex = 1:N_figures
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex = 1:sp1*sp2
%         if dex1 <= N_outputs
%             subplot(sp1,sp2,subdex)
%             % Plot data points
%             plot(exp(Y(:,dex1)),exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color_1);
%             xlabel(['Actual ', output_names{dex1}])
%             ylabel(['Predicted ', output_names{dex1}])
%             title(['R^2 = ',num2str(R2each(dex1),4)])
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             % Plot identity line
%             ylim_ind = get(gca,'ylim') ;
%             xlim_ind = get(gca,'xlim') ;
%             minpoint = min([ylim_ind(1),xlim_ind(1)]);
%             maxpoint = max([ylim_ind(2),xlim_ind(2)]);
%             hold on
%             plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
%             dex1 = dex1+1;
%         end
%     end
% end

dex2 = 1;
for figdex2 = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:sp1*sp2
        if dex2 <= N_outputs
            subplot(sp1,sp2,subdex2)
            bar(Bpls(:,dex2),'FaceColor',color_1)
            title(output_names(dex2))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex2 = dex2 + 1;
        end
    end
end

Bpls_1 = Bpls;

disp('Mouse')
good_count
avg_R2
disp('*****************')

%% RABBIT

all_outputs = all_outputs_2(:,output_selection);

% Outputs - population level
%all_outputs_mean = mean(all_outputs);
%all_outputs_std_dev = std(all_outputs);

% Check basic properties, and define X and Y
% Exclude cells with outputs = 0
%good_trials = all_outputs(:,1)>0;
    good_trials = (all_outputs(:,1)>0) .* (all_outputs(:,2)<400); % APD<400ms
good_trials_log = logical(good_trials);
% Good count: number of good trials
good_count = sum(good_trials_log);
% Good_parameters: array with parameters from good trials only
good_parameters = all_parameters(good_trials_log,:);
X = log(good_parameters);
% Good_outputs: array with parameters from good trials only
good_outputs = all_outputs(good_trials_log,:);
Y = log(good_outputs);

%all_outputs_mean = mean(good_outputs);
%all_outputs_std_dev = std(good_outputs);

% Istograms
dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_2)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% % PLS - svds algorithm (2010)
% [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,Yjack,R2x,R2y,...
%           RESSy,PRESSy,Q2,r2y_random,rv_random,Yhat4Press,Yhat4Ress] = ...
%           PLS_jack_svds(X,Y,rank(X));

% N_pars, N_outputs, N_trials, good_count

% Calculate agreement of values predicted by regression (Yhat = Bpls*X) 
% with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

avg_R2 = mean(R2each);

% R2
figure(100),set(gcf,'color','w','Position',[50,100,1500,750])
subplot(3,1,2)
bar(R2each,'FaceColor',color_2)
set(gca,'box','off','tickdir','out','fontsize',12)
set(gca,'XTick',1:N_outputs)
set(gca,'XTickLabel',output_names)
set(gca,'XLim',[0 N_outputs+1])
ylim([0 1])
rotateXLabels( gca(), 90)
%title('R^2 values')
title(['R^2 values (mean = ',num2str(avg_R2,3),')'])

% dex1 = 1;
% for figdex = 1:N_figures
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex = 1:sp1*sp2
%         if dex1 <= N_outputs
%             subplot(sp1,sp2,subdex)
%             % Plot data points
%             plot(exp(Y(:,dex1)),exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color_2);
%             xlabel(['Actual ', output_names{dex1}])
%             ylabel(['Predicted ', output_names{dex1}])
%             title(['R^2 = ',num2str(R2each(dex1),4)])
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             % Plot identity line
%             ylim_ind = get(gca,'ylim') ;
%             xlim_ind = get(gca,'xlim') ;
%             minpoint = min([ylim_ind(1),xlim_ind(1)]);
%             maxpoint = max([ylim_ind(2),xlim_ind(2)]);
%             hold on
%             plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
%             dex1 = dex1+1;
%         end
%     end
% end

dex2 = 1;
for figdex2 = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:sp1*sp2
        if dex2 <= N_outputs
            subplot(sp1,sp2,subdex2)
            bar(Bpls(:,dex2),'FaceColor',color_2)
            title(output_names(dex2))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex2 = dex2 + 1;
        end
    end
end

Bpls_2 = Bpls;

disp('Rabbit')
good_count
avg_R2
disp('*****************')

%% HUMAN

all_outputs = all_outputs_3(:,output_selection);

% Outputs - population level
%all_outputs_mean = mean(all_outputs);
%all_outputs_std_dev = std(all_outputs);

% Check basic properties, and define X and Y
% Exclude cells with outputs = 0
%good_trials = all_outputs(:,1)>0;
    good_trials = (all_outputs(:,1)>0) .* (all_outputs(:,2)<400); % APD<400ms
good_trials_log = logical(good_trials);
% Good count: number of good trials
good_count = sum(good_trials_log);
% Good_parameters: array with parameters from good trials only
good_parameters = all_parameters(good_trials_log,:);
X = log(good_parameters);
% Good_outputs: array with parameters from good trials only
good_outputs = all_outputs(good_trials_log,:);
Y = log(good_outputs);

%all_outputs_mean = mean(good_outputs);
%all_outputs_std_dev = std(good_outputs);

% Istograms
dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_3)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% % PLS - svds algorithm (2010)
% [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,Yjack,R2x,R2y,...
%           RESSy,PRESSy,Q2,r2y_random,rv_random,Yhat4Press,Yhat4Ress] = ...
%           PLS_jack_svds(X,Y,rank(X));

% N_pars, N_outputs, N_trials, good_count

% Calculate agreement of values predicted by regression (Yhat = Bpls*X) 
% with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

avg_R2 = mean(R2each);

% R2
figure(100),set(gcf,'color','w','Position',[50,100,1500,750])
subplot(3,1,3)
bar(R2each,'FaceColor',color_3)
set(gca,'box','off','tickdir','out','fontsize',12)
set(gca,'XTick',1:N_outputs)
set(gca,'XTickLabel',output_names)
set(gca,'XLim',[0 N_outputs+1])
ylim([0 1])
rotateXLabels( gca(), 90)
%title('R^2 values')
title(['R^2 values (mean = ',num2str(avg_R2,3),')'])

% dex1 = 1;
% for figdex = 1:N_figures
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex = 1:sp1*sp2
%         if dex1 <= N_outputs
%             subplot(sp1,sp2,subdex)
%             % Plot data points
%             plot(exp(Y(:,dex1)),exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color_3);
%             xlabel(['Actual ', output_names{dex1}])
%             ylabel(['Predicted ', output_names{dex1}])
%             title(['R^2 = ',num2str(R2each(dex1),4)])
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             % Plot identity line
%             ylim_ind = get(gca,'ylim') ;
%             xlim_ind = get(gca,'xlim') ;
%             minpoint = min([ylim_ind(1),xlim_ind(1)]);
%             maxpoint = max([ylim_ind(2),xlim_ind(2)]);
%             hold on
%             plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
%             dex1 = dex1+1;
%         end
%     end
% end

dex2 = 1;
for figdex2 = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:sp1*sp2
        if dex2 <= N_outputs
            subplot(sp1,sp2,subdex2)
            bar(Bpls(:,dex2),'FaceColor',color_3)
            title(output_names(dex2))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex2 = dex2 + 1;
        end
    end
end

Bpls_3 = Bpls;

disp('Human')
good_count
avg_R2
disp('*****************')

%%
sp1 = 2; sp2 = 1;
max_panels = 2;
N_figures = ceil(N_outputs/max_panels);

dex3 = 1;
for figdex2 = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:sp1*sp2
        if dex3 <= N_outputs
            subplot(sp1,sp2,subdex2)
            
            Bpls = [Bpls_1(:,dex3) Bpls_2(:,dex3) Bpls_3(:,dex3)];
            hb = bar(Bpls);%,'FaceColor',color)
            set(hb(1), 'FaceColor',color_1)%, 'EdgeColor',color_nSR_1_e)
            set(hb(2), 'FaceColor',color_2)%, 'EdgeColor',color_cAF_1_e)
            set(hb(3), 'FaceColor',color_3)%, 'EdgeColor',color_cAF_1_e)
            title(output_names(dex3))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex3 = dex3 + 1;
        end
    end
end
