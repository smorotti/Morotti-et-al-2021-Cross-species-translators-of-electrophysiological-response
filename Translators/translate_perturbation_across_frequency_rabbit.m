%% Cross-pacing rate translation of the effect of drugs administration
% (ion channel blockers or ISO) observed in rabbit myocytes

close all
clear
clc

%% Perturbation selection
drug_index = 1; % 1 INaL, 2 Ito, 3 Ik1, 4 Ikr, 5 Iks, 6 Ikca, 7 ISO

% List of perturbations:
% 1) INaL inhibition (GS-967 pretreatment) 1 uM
% 2) Ito inhibition (4-AP pretreatment) 5 mM
% 3) IK1 inhibition (PA-6 pretreatment) 200 nM
% 4) IKr inhibition (E-4031 pretreatment) 1 uM
% 5) IKs inhibition (HMR-1556 pretreatment) 1 uM
% 6) IKCa inhibition (Apamin pretreatment) 100 nM
% 7) ISO administration (30 nM)

%% Load and plot experimental data (Hegyi et al. 2021)
disp('Rabbit - Exp data from Hegyi et al. 2021')

color_rabbit_1Hz = [0.65 0.65 0.65];
color_rabbit = [0.85 0.33 0.1];

Freq = [0.5 1 2 3];

% Control
APD50_data = [163.2582754	9.284659681	18;	170.1345514	7.371768043	18;	151.2051578	9.072109889	18;	125.3984221	6.361444887	18];
APD90_data = [211.7963255	9.763974846	18;	221.5581439	8.001801947	18;	206.7929628	9.97623836	18;	182.9544444	6.826626298	18];

APD50_control = APD50_data(:,1); APD50_SEM_control = APD50_data(:,2); APD50_STD_control = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_control = APD90_data(:,1); APD90_SEM_control = APD90_data(:,2); APD90_STD_control = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_control = 18;

% INaL inhibition (GS-967 pretreatment)
APD50_data = [118.241086	16.58308567	4;	125.5849715	8.266857887	7;	129.7720005	13.12310024 4;	111.113333	9.221632951	4];
APD90_data = [160.9841887	11.42294522	4;	164.0821719	7.121177946	7;	166.5623997	9.344941888	4;	148.2331996	5.445881986	4];

APD50_inal = APD50_data(:,1); APD50_SEM_inal = APD50_data(:,2); APD50_STD_inal = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_inal = APD90_data(:,1); APD90_SEM_inal = APD90_data(:,2); APD90_STD_inal = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_inal = 4;

% Ito inhibition (4-AP pretreatment)
APD50_data = [230.9499463	10.21168587	8;	199.3782857	6.602695984	16;	172.0098667	8.95065952	6;	135.9499338	7.910797747	6];
APD90_data = [275.0768542	10.3110336	8;	247.7622858	6.788192816	16;	215.0495336	9.873405458	6;	178.1584667	8.299868032	6];

APD50_ito = APD50_data(:,1); APD50_SEM_ito = APD50_data(:,2); APD50_STD_ito = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_ito = APD90_data(:,1); APD90_SEM_ito = APD90_data(:,2); APD90_STD_ito = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_ito = 6;%8;

% IK1 inhibition (PA-6 pretreatment)
APD50_data = [174.8469471	17.69209532	11;	179.798764	14.34862108	11;	165.6217816	10.95720211	11;	138.369818	8.031444782	11];
APD90_data = [255.3053426	16.57819902	11;	257.21142	12.92680855	11;	231.82782	10.68066117	11;	196.38782	7.70120974	11];

APD50_ik1 = APD50_data(:,1); APD50_SEM_ik1 = APD50_data(:,2); APD50_STD_ik1 = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_ik1 = APD90_data(:,1); APD90_SEM_ik1 = APD90_data(:,2); APD90_STD_ik1 = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_ik1 = 10;

% IKr inhibition (E-4031 pretreatment)
APD50_data = [209.8423873	7.839285195	9;	203.49804	6.71191419	10;	180.4156802	6.713989007	10;	141.3352003	6.181357268	10];
APD90_data = [266.8963775	5.839530165	9;	260.10124	3.748317632	10;	236.2891198	5.80666467	10;	193.7738798	3.818453697	10];

APD50_ikr = APD50_data(:,1); APD50_SEM_ikr = APD50_data(:,2); APD50_STD_ikr = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_ikr = APD90_data(:,1); APD90_SEM_ikr = APD90_data(:,2); APD90_STD_ikr = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_ikr = 9;

% IKs inhibition (HMR-1556 pretreatment)
APD50_data = [159.2978302	9.786793325	10;	166.27536	9.880192926	10;	163.5300799	5.817173784	10;	142.6277602	5.576759108	10];
APD90_data = [209.76083	10.95215207	10;	215.05304	10.33047524	10;	209.61532	5.703759419	10;	186.0128	5.072987137	10];

APD50_iks = APD50_data(:,1); APD50_SEM_iks = APD50_data(:,2); APD50_STD_iks = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_iks = APD90_data(:,1); APD90_SEM_iks = APD90_data(:,2); APD90_STD_iks = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_iks = 10;

% IKCa inhibition (Apamin pretreatment)
APD50_data = [166.7460084	10.3661082	7;	172.6949999	7.216544081	8;	153.7113001	6.619783622	8;	123.0636505	5.622938692	8];
APD90_data = [216.00018	9.745737495	7;	223.2250002	5.832823164	8;	202.3856	4.982563417	8;	170.29495	3.824706403	8];

APD50_ikca = APD50_data(:,1); APD50_SEM_ikca = APD50_data(:,2); APD50_STD_ikca = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_ikca = APD90_data(:,1); APD90_SEM_ikca = APD90_data(:,2); APD90_STD_ikca = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_ikca = 7;

% ISO
APD50_data = [146.7841487	10.66501524	6;	150.5785142	6.6650293	7;	137.7037999	5.855803388	10;	115.4536667	3.989066413	6];
APD90_data = [186.43485	8.856615372	6;	188.19817	5.262625723	7;	173.52693	3.916940618	10;	147.9958	2.349790371	6];

APD50_iso = APD50_data(:,1); APD50_SEM_iso = APD50_data(:,2); APD50_STD_iso = APD50_data(:,2).*sqrt(APD50_data(:,3));
APD90_iso = APD90_data(:,1); APD90_SEM_iso = APD90_data(:,2); APD90_STD_iso = APD90_data(:,2).*sqrt(APD90_data(:,3));
n_iso = 6;

shift_ctrl_90 = -0.025;
shift_ctrl_50 = +0.025;
shift_drug_90 = 0;
shift_drug_50 = +0.05;

if drug_index == 1
    disp('Cross-frequency translation of INaL block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_inal,APD90_STD_inal,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_inal,APD50_STD_inal,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 300])
    title('Rabbit - INaL block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 2
    disp('Cross-frequency translation of Ito block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_ito,APD90_STD_ito,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_ito,APD50_STD_ito,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 350])
    title('Rabbit - Ito block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 3
    disp('Cross-frequency translation of IK1 block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_ik1,APD90_STD_ik1,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_ik1,APD50_STD_ik1,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 350])
    title('Rabbit - IK1 block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 4
    disp('Cross-frequency translation of IKr block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_ikr,APD90_STD_ikr,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_ikr,APD50_STD_ikr,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 350])
    title('Rabbit - IKr block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 5
    disp('Cross-frequency translation of IKs block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_iks,APD90_STD_iks,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_iks,APD50_STD_iks,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 350])
    title('Rabbit - IKs block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 6
    disp('Cross-frequency translation of IK,Ca block')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_ikca,APD90_STD_ikca,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_ikca,APD50_STD_ikca,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'),xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 300]) % 300 supplement
    title('Rabbit - IK,Ca block')
    legend('APD90 ctrl','APD50 ctrl','APD90 drug','APD50 drug')
elseif drug_index == 7
    disp('Cross-frequency translation of ISO administration')
    
    figure,set(gcf,'color','w'),hold on
    errorbar(Freq+shift_ctrl_90,APD90_control,APD90_STD_control,'k--','LineWidth',1)
    errorbar(Freq+shift_ctrl_50,APD50_control,APD50_STD_control,'r--','LineWidth',1)
    errorbar(Freq+shift_drug_90,APD90_iso,APD90_STD_iso,'Color',[0.5 0.5 0.5],'LineWidth',1)
    errorbar(Freq+shift_drug_50,APD50_iso,APD50_STD_iso,'Color',[0.850980392156863 0.329411764705882 0.101960784313725],'LineWidth',1)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Frequency (Hz)'), xlim([0.25 3.25])
    ylabel('APD (ms)'),ylim([0 300])
    title('Rabbit - ISO')
    legend('APD90 ctrl','APD50 ctrl','APD90 iso','APD50 iso')
end

%% Translation
switch drug_index
    case 1
        APD50_drug = APD50_inal; APD50_SEM_drug = APD50_SEM_inal;
        APD90_drug = APD90_inal; APD90_SEM_drug = APD90_SEM_inal;
        n_drug = n_inal;
    case 2
        APD50_drug = APD50_ito; APD50_SEM_drug = APD50_SEM_ito;
        APD90_drug = APD90_ito; APD90_SEM_drug = APD90_SEM_ito;
        n_drug = n_ito;
    case 3
        APD50_drug = APD50_ik1; APD50_SEM_drug = APD50_SEM_ik1;
        APD90_drug = APD90_ik1; APD90_SEM_drug = APD90_SEM_ik1;
        n_drug = n_ik1;
    case 4
        APD50_drug = APD50_ikr; APD50_SEM_drug = APD50_SEM_ikr;
        APD90_drug = APD90_ikr; APD90_SEM_drug = APD90_SEM_ikr;
        n_drug = n_ikr;
    case 5
        APD50_drug = APD50_iks; APD50_SEM_drug = APD50_SEM_iks;
        APD90_drug = APD90_iks; APD90_SEM_drug = APD90_SEM_iks;
        n_drug = n_iks;
    case 6
        APD50_drug = APD50_ikca; APD50_SEM_drug = APD50_SEM_ikca;
        APD90_drug = APD90_ikca; APD90_SEM_drug = APD90_SEM_ikca;
        n_drug = n_ikca;
    case 7
        APD50_drug = APD50_iso; APD50_SEM_drug = APD50_SEM_iso;
        APD90_drug = APD90_iso; APD90_SEM_drug = APD90_SEM_iso;
        n_drug = n_iso;
end

pos_freq_input = 2; % 1 Hz
pos_freq_output = [1 3 4]; % 0.5, 2, 3 Hz

n_pop_fitting = zeros(1,length(pos_freq_output));
avg_R2_fitting = zeros(1,length(pos_freq_output));
avg_R2_validation = zeros(1,length(pos_freq_output));

for ii = 1:length(pos_freq_output)
    % Load input data
    load outputs_matrix_1500_0p1_300s_rabbit_COMMON % 1 Hz
    all_outputs_rabbit_control = all_outputs;
    
    switch pos_freq_output(ii)
        case 1
            disp('Control-1 to Control-0p5')
            
            load outputs_matrix_1500_0p1_300s_rabbit_COMMON_0p5Hz
            all_outputs_rabbit_control_0p5Hz = all_outputs;
            
            Xblock = all_outputs_rabbit_control;
            Yblock = all_outputs_rabbit_control_0p5Hz;
        case 3
            disp('Control-1 to Control-2')
            
            load outputs_matrix_1500_0p1_300s_rabbit_COMMON_2Hz
            all_outputs_rabbit_control_2Hz = all_outputs;
            
            Xblock = all_outputs_rabbit_control;
            Yblock = all_outputs_rabbit_control_2Hz;
        case 4
            disp('Control-1 to Control-3')
            
            load outputs_matrix_1500_0p1_300s_rabbit_COMMON_3Hz
            all_outputs_rabbit_control_3Hz = all_outputs;
            
            Xblock = all_outputs_rabbit_control;
            Yblock = all_outputs_rabbit_control_3Hz;
    end
    
    % Select outputs
    output_selection = [5 7]; % 2 outputs: APD90, APD50
    
    Xblock = Xblock(:,output_selection);
    Yblock = Yblock(:,output_selection);
    
    output_names = output_names(output_selection);
    output_units = output_units(output_selection);
    
    % N_outputs = length(output_names);
    
    output_names_X = output_names;
    output_units_X = output_units;
    
    N_outputs_X = length(output_names_X);
    
    output_names_Y = output_names;
    output_units_Y = output_units;
    
    % Check basic properties, and define X and Y
    % Exclude cells with outputs = 0
    good_trials = (Xblock(:,1)>0).*(Yblock(:,1)>0);
    good_trials_log = logical(good_trials);
    % Good count: number of good trials
    good_count_total = sum(good_trials_log);
    % Good_outputs: array with parameters from good trials only
    good_outputs_X = Xblock(good_trials_log,:);
    good_outputs_Y = Yblock(good_trials_log,:);
    %X = log(good_outputs_X);
    Y = log(good_outputs_Y);
    
    % % No check
    % good_count = N_trials;
    % X = log(Xblock);
    % Y = log(Yblock);
    
    %[N_trials N_outputs_X] = size(X);
    [N_trials N_outputs_Y] = size(Y);
    
    % Fitting vs Validation Groups
    Nval = 400;
    
    actual_inputs = good_outputs_X(end-Nval+1:end,:);
    actual_outputs = good_outputs_Y(end-Nval+1:end,:);
    
    good_count = good_count_total-Nval;
    X = log(good_outputs_X(1:end-Nval,:));
    Y = log(good_outputs_Y(1:end-Nval,:));
    
    % Plotting options
    color = [0 0 0]; % black
    
    plot_matrix = 0;
    plot_fitting = 0;
    plot_validation = 0;
    
    % Call the PLS routine
    N_fitting = good_count;
    
    % PLS - nipals algorithm (2003)
    [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
        PLS_nipals(X,Y,rank(X));
    % N_pars, N_outputs, N_trials, good_count
    
    % Goodness of fit - R^2
    % Fraction of the variance in the dependent variable which is explained by the model
    
    % Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
    % Assessment on log-transformed values
    SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
    SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
    R2each = SSYR./SSYT;
    avg_R2_fit = mean(R2each);
    
    % Assessment on (normal) values
    R2ord_fit = zeros(1,N_outputs_Y);
    R2adj_fit = zeros(1,N_outputs_Y);
    rmse_fit = zeros(1,N_outputs_Y);
    for i = 1:N_outputs_Y
        mdl = fitlm(exp(Y(:,i)),exp(Yhat(:,i)));
        R2ord_fit(i) = mdl.Rsquared.Ordinary;
        R2adj_fit(i) = mdl.Rsquared.Adjusted;
        rmse_fit(i) = mdl.RMSE;
    end
    R2ord_fit;
    R2adj_fit;
    
    R2_fit = R2adj_fit; % Values plotted in figures
    avg_R2calc_fit = mean(R2_fit);
    
    % Residual Standard Deviation
    % Standard deviation for a normal distribution, centered on the predicted regression line,
    % representing the distribution of actually observed values
    oSD = zeros(1,N_outputs_Y);
    rSD = zeros(1,N_outputs_Y);
    for i = 1:N_outputs_Y
        oSD(i) = std(exp(Y(:,i)));
        rSD(i) = sqrt(sum((exp(Yhat(:,i)) - exp(Y(:,i)) ).^2) / (good_count-2));
    end
    
    % Plot regression coefficients
    if plot_matrix == 1
        figure; set(gcf,'color','w')
        imagesc(Bpls'); colormap jet;
        set(gca,'box','off','tickdir','out','fontsize',14)
        set(gca,'YDir','normal')
        title('Regression coefficients');
        xlabel('Outputs 1');
        ylabel('Outputs 2');
        set(gca,'XTick',(1:N_outputs_X))
        set(gca,'XTickLabel',output_names_X)
        set(gca,'YTickLabel',output_names_Y)
        set(gca,'YTick',(1:N_outputs_Y))
        rotateXLabels(gca(), 90)
        colorbar
    end
    
    % Scatter plot
    N_figures = ceil(N_outputs_Y/10);
    
    if plot_fitting == 1
        dex1 = 1;
        for figdex = 1:N_figures
            figure
            set(gcf,'color','w','Position',[50,100,1500,750])
            for i = 1:10
                if dex1 <= N_outputs_Y
                    subplot(2,5,i),hold on
                    % Plot data points
                    factor = 1;
                    if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                    if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                    plot(factor*exp(Y(:,dex1)),factor*exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color)
                    xlabel(['Actual ',output_names_Y{dex1}])
                    ylabel(['Predicted ',output_names_Y{dex1}])
                    title(['R^2 = ',num2str(R2_fit(dex1),4)])
                    set(gca,'box','off','tickdir','out','fontsize',12)
                    % Plot identity line
                    ylim_ind = get(gca,'ylim') ;
                    xlim_ind = get(gca,'xlim') ;
                    minpoint = min([ylim_ind(1),xlim_ind(1)]);
                    maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                    plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                    xlim([minpoint, maxpoint])
                    ylim([minpoint, maxpoint])
                    dex1 = dex1+1;
                end
            end
        end
    end
    
    % Validation
    N_validation = Nval;
    
    predicted_outputs = actual_outputs*0;
    for i = 1:Nval
        cell_index = i;
        inputs = actual_inputs(cell_index,:);
        x = log(inputs);
        xz = (x-mean(X))./std(X);
        
        yz = xz*Bpls;
        %    yz = xz;
        y = yz.*std(Y)+mean(Y);
        predicted_outputs(i,:) = exp(y);
    end
    
    R2ord = zeros(1,N_outputs_Y);
    R2adj = zeros(1,N_outputs_Y);
    rmse_val = zeros(1,N_outputs_Y);
    for i = 1:N_outputs_Y
        mdl = fitlm(actual_outputs(:,i),predicted_outputs(:,i));
        R2ord(i) = mdl.Rsquared.Ordinary;
        R2adj(i) = mdl.Rsquared.Adjusted;
        rmse_val(i) = mdl.RMSE;
    end
    R2ord;
    R2adj;
    R2 = R2adj;
    avg_R2_val = mean(R2);
    
    % Residual Standard Deviation
    % Standard deviation for a normal distribution, centered on the predicted regression line,
    % representing the distribution of actually observed values
    oSD_val = zeros(1,N_outputs_Y);
    rSD_val = zeros(1,N_outputs_Y);
    for i = 1:N_outputs_Y
        oSD_val(i) = std(actual_outputs(:,i));
        rSD_val(i) = sqrt(sum((predicted_outputs(:,i) - actual_outputs(:,i) ).^2) / (N_validation-2));
    end
    
    %plot_validation = 1;
    if plot_validation == 1
        dex1 = 1;
        for figdex = 1:N_figures
            figure
            set(gcf,'color','w','Position',[50,100,1500,750])
            for i = 1:10
                if dex1 <= N_outputs_Y
                    subplot(2,5,i),hold on
                    % Plot data points
                    factor = 1;
                    if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                    if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                    plot(factor*actual_outputs(:,dex1),factor*predicted_outputs(:,dex1),'Marker','o','LineStyle','none')%,'Color',color)
                    xlabel(['Actual ',output_names_Y{dex1}])
                    ylabel(['Predicted ',output_names_Y{dex1}])
                    title(['R^2 = ',num2str(R2(dex1),4)])
                    set(gca,'box','off','tickdir','out','fontsize',14)
                    % Plot identity line
                    ylim_ind = get(gca,'ylim') ;
                    xlim_ind = get(gca,'xlim') ;
                    minpoint = min([ylim_ind(1),xlim_ind(1)]);
                    maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                    plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                    xlim([minpoint, maxpoint])
                    ylim([minpoint, maxpoint])
                    dex1 = dex1+1;
                end
            end
        end
    end
    
    n_pop_fitting(ii) = N_fitting;
    avg_R2_fitting(ii) = avg_R2_fit;
    avg_R2_validation(ii) = avg_R2_val;
    
    % Control
    rabbit_control = [APD90_control(pos_freq_input) APD50_control(pos_freq_input)]; % control rabbit 1 Hz
    rabbit_control_sem = [APD90_SEM_control(pos_freq_input) APD50_SEM_control(pos_freq_input)];
    rabbit_control_std = rabbit_control_sem*sqrt(n_control);
    rabbit_control_xHz = [APD90_control(pos_freq_output(ii)) APD50_control(pos_freq_output(ii))]; % control rabbit 0.5, 2, 3 Hz
    rabbit_control_xHz_sem = [APD90_SEM_control(pos_freq_output(ii)) APD50_SEM_control(pos_freq_output(ii))];
    rabbit_control_xHz_std = rabbit_control_xHz_sem*sqrt(n_control);
    % Drug
    rabbit_drug = [APD90_drug(pos_freq_input) APD50_drug(pos_freq_input)]; % drug rabbit 1 Hz
    rabbit_drug_sem = [APD90_SEM_drug(pos_freq_input) APD50_SEM_drug(pos_freq_input)];
    rabbit_drug_std = rabbit_drug_sem*sqrt(n_drug);
    rabbit_drug_xHz = [APD90_drug(pos_freq_output(ii)) APD50_drug(pos_freq_output(ii))]; % drug rabbit 0.5, 2, 3 Hz
    rabbit_drug_xHz_sem = [APD90_SEM_drug(pos_freq_output(ii)) APD50_SEM_drug(pos_freq_output(ii))];
    rabbit_drug_xHz_std = rabbit_drug_xHz_sem*sqrt(n_drug);
    
    % Ratio (drug/drug-free)
    ratio_rabbit_std = (rabbit_drug./rabbit_control).*sqrt((rabbit_drug_std./rabbit_drug).^2+(rabbit_control_std./rabbit_control).^2);
    ratio_rabbit_xHz_std = (rabbit_drug_xHz./rabbit_control_xHz).*sqrt((rabbit_drug_xHz_std./rabbit_drug_xHz).^2+(rabbit_control_xHz_std./rabbit_control_xHz).^2);
    ratio_rabbit_sem = ratio_rabbit_std./sqrt(n_drug);
    ratio_rabbit_xHz_sem = ratio_rabbit_xHz_std./sqrt(n_drug);
    
    % Simulated data (control)
    rabbit_sim = mean(good_outputs_X);
    rabbit_sim_std = std(good_outputs_X);
    rabbit_xHz_sim = mean(good_outputs_Y);
    rabbt_xHz_sim_std = std(good_outputs_Y);
    
    % Predicition (with normalized drug-effect)
    actual_outputs = rabbit_drug_xHz./rabbit_control_xHz.*(mean(good_outputs_Y));
    actual_inputs = rabbit_drug./rabbit_control.*(mean(good_outputs_X));
    rabbit_drug_xHz_std_scaled = ratio_rabbit_xHz_std.*(mean(good_outputs_Y));
    rabbit_drug_std_scaled = ratio_rabbit_std.*(mean(good_outputs_X));
    rabbit_drug_sem_scaled = rabbit_drug_xHz_std_scaled./sqrt(n_drug);
    rabbit_drug_xHz_sem_scaled = rabbit_drug_std_scaled./sqrt(n_drug);
    
    inputs = actual_inputs;
    x = log(inputs);
    xz = (x-mean(X))./std(X);
    
    yz = xz*Bpls;
    y = yz.*std(Y)+mean(Y);
    predicted_outputs = exp(y);
    
    figure(100),set(gcf,'color','w'),hold on
    subplot(1,3,ii), hold on
    b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
    set(b(1),'FaceColor',color_rabbit_1Hz); set(b(3),'FaceColor',color_rabbit);
    errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[actual_inputs;predicted_outputs;actual_outputs]',[rabbit_drug_std_scaled; 0 0 ;rabbit_drug_xHz_std_scaled]','k.')
    set(gca,'box','off','tickdir','out','fontsize',14)
    %legend('Rabbit 1-Hz','1-to-x Hz','Rabbit x-Hz')
    set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
end

figure(100)
subplot(1,3,1),title('0.5 Hz')
ylabel('"Scaled" APD values after drug (ms)')
subplot(1,3,2),title('2 Hz')
subplot(1,3,3),title('3 Hz')
legend('1-Hz data','Predicted','Actual')
