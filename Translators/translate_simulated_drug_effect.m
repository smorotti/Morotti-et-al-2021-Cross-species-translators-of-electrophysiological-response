%% Mouse-to-Human and Rabbit-to-Human translation of data obtained
% with the baseline models in control condition or upon selective (50%)
% reduction of GNaL, Gtof, GKr, GK1, GCaL, vNCX, vSERCA, or vRel.

close all
clear
clc

%% Selection of the condition to test
% Change current_index_SA to select the desired condition (see list below)
current_index_SA = 17;

% Possible options for current_index_SA:
% 0) Control (drug-free)
% 2) 50% GNaL
% 5) 50% Gtof
% 8) 50% GKr
% 13) 50% GK1
% 17) 50% GCaL
% 19) 50% vNCX
% 21) 50% vSERCA
% 22) 50% vRel

%% Parameters
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak
parameter_names = {'GNa' 'GNaL' 'GNaB' 'vNKA' 'Gtof'...
    'Gtos' 'GKs' 'GKr' 'GKur1' 'GKur2'...
    'Gss' 'GKp' 'GK1' 'GCFTR' 'GClCa'...
    'GClB' 'GCaL' 'GCaB' 'vNCX' 'vPMCA'...
    'vSERCA' 'vRel' 'vLeak'} ;

% Currents/Transporters to block
block_index = [2 5 8 13 17 19 21 22];

if current_index_SA == 0
    disp('Drug-free')
else
    current_index = find(block_index == current_index_SA);
    if isempty(current_index) == 1
        disp('Wrong selection!')
        exit_execution % fake command
    else
        disp('Parameter selected:')
        disp(parameter_names{current_index_SA})
    end
end

%% Loading simulated dataset (outputs)
load outputs_current_block_human % outputs_human_block
load outputs_current_block_rabbit % outputs_rabbit_block
load outputs_current_block_mouse % outputs_mouse_block

output_names_original = {'UV', 'APpeak', '|MDP|', 'APamp', 'APD90',...
    'APD70', 'APD50', 'APD30', 'CaTmax', 'CaTmin',...
    'CaTamp', 'CaTttp', 'CaTt50', 'CaTtau', 'Na',...
    'CaSRmax', 'CaSRmin', 'CaSRamp'};

% output_units = {'mV/ms', 'mV', 'mV', 'mV', 'ms',...
%     'ms', 'ms', 'ms', 'mM', 'mM',...
%     'mM', 'ms', 'ms', 'ms', 'mM',...
%     'mM', 'mM', 'mM'};

%% Cycle
% Colors
color_human = [0 0.45 0.74];
color_rabbit = [0.85 0.33 0.1];
color_mouse = [0.47 0.67 0.19];

output_index_array = [1 2 3 4];

all_R2_output_mouse = zeros(length(output_index_array),10);
all_R2_output_rabbit = zeros(length(output_index_array),10);

all_R2_output_index_mouse = 1;
all_R2_output_index_rabbit = 2;

all_predictions_from_mouse = zeros(length(output_index_array),10);
all_predictions_from_rabbit = zeros(length(output_index_array),10);

for ii = 1:length(output_index_array)
    output_index = output_index_array(ii);
    
    % Output selection
    % 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
    % 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
    % 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
    % 16) CaSRmax 17) CaSRmin 18) CaSRamp
    
    if output_index == 0 % ALL
        output_selection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    % ALL OUTPUTS
    elseif output_index == 1
        output_selection = [1 3 4 5 7 10 11 12 13 14]; % 10
    % UV	MDP	APamp	APD90	APD50	CaTmin	CaTamp	CaTttp	CaTt50	CaTtau
    elseif output_index == 2
        output_selection = [5 7 10 11 13 14]; % 6
    % APD90	APD50	CaTmin	CaTamp	CaTt50	CaTtau
    elseif output_index == 3
        output_selection = [5 7 13 14]; % 4
    % APD90	APD50	CaTt50	CaTtau
    elseif output_index == 4
        output_selection = [5 7]; % 2
    % APD90	APD50
    end
    
    X = ['# features ',num2str(length(output_selection))];
    disp(X)

    % Populations selection: Mouse-to-Human, Rabbit-to-Human
    pop_index = [all_R2_output_index_mouse all_R2_output_index_rabbit];

    n_pop_fitting = zeros(1,length(pop_index));
    avg_R2_fitting = zeros(1,length(pop_index));
    avg_R2_validation = zeros(1,length(pop_index));  

    figure(1),set(gcf,'color','w','Position',[50,100,1500,750])
    for j = 1:length(pop_index)
        %j
        translator_index = pop_index(j);
        
        % Loading population-based simulated data
        load outputs_matrix_1500_0p1_300s_human_COMMON
        all_outputs_human_control = all_outputs;

        load outputs_matrix_1500_0p1_300s_rabbit_COMMON
        all_outputs_rabbit_control = all_outputs;

        load outputs_matrix_1500_0p1_300s_mouse_COMMON
        all_outputs_mouse_control = all_outputs;
        
        switch translator_index
            case 1
                %disp('Mouse-to-Human')
                Xblock = all_outputs_mouse_control;
                Yblock = all_outputs_human_control;
            case 2
                %disp('Rabbit-to-Human')
                Xblock = all_outputs_rabbit_control;
                Yblock = all_outputs_human_control;
        end

        % Select outputs
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

        % with R^2 = 0.75, the model explains approximately 75% of the variability in the predicted variable
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

        n_pop_fitting(j) = N_fitting;
        avg_R2_fitting(j) = avg_R2_fit;
        avg_R2_validation(j) = avg_R2_val;

        if translator_index == all_R2_output_index_mouse
            all_R2_output_mouse(ii,1:length(R2)) = R2;
        end
        if translator_index == all_R2_output_index_rabbit
            all_R2_output_rabbit(ii,1:length(R2)) = R2;
        end
        
        % Application to simulated data
        if current_index_SA == 0 % Drug-free
            if translator_index == all_R2_output_index_mouse
                actual_inputs = outputs_mouse_baseline(output_selection);
            end
            if translator_index == all_R2_output_index_rabbit
                actual_inputs = outputs_rabbit_baseline(output_selection);
            end
            actual_outputs = outputs_human_baseline(output_selection);
        else
            if translator_index == all_R2_output_index_mouse
                actual_inputs = outputs_mouse_block(current_index,output_selection);
            end
            if translator_index == all_R2_output_index_rabbit
                actual_inputs = outputs_rabbit_block(current_index,output_selection);
            end
            actual_outputs = outputs_human_block(current_index,output_selection);
        end
        
        output_names_selected = output_names_original(output_selection);
        
        inputs = actual_inputs;
        x = log(inputs);
        xz = (x-mean(X))./std(X);

        yz = xz*Bpls;
        y = yz.*std(Y)+mean(Y);
        predicted_outputs = exp(y);
        
        position_index_0 = 0; color_input = color_mouse;
        if translator_index == all_R2_output_index_rabbit
            position_index_0 = length(output_index_array);
            color_input = color_rabbit;
        end
        figure(1)%,set(gcf,'color','w','Position',[50,100,1500,750])
        subplot(2,length(output_index_array),position_index_0+ii),hold on
        b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
        set(b(1),'FaceColor',color_input); set(b(3),'FaceColor',color_human);
        set(gca,'box','off','tickdir','out','fontsize',8)
        xlabel('Outputs'),ylabel('Values')
        xticks(1:length(actual_outputs));
        xticklabels(output_names_selected)
        
        if translator_index == all_R2_output_index_mouse
            all_predictions_from_mouse(ii,1:length(output_selection)) = predicted_outputs;
        end
        if translator_index == all_R2_output_index_rabbit
            all_predictions_from_rabbit(ii,1:length(output_selection)) = predicted_outputs;
        end
    end
end

figure(1)
if current_index_SA == 0
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, control')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, control')
elseif current_index_SA == 2
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% INaL')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% INaL')
elseif current_index_SA == 5
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% Itof')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% Itof')
elseif current_index_SA == 8
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% IKr')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% IKr')
elseif current_index_SA == 13
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% IK1')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% IK1')
elseif current_index_SA == 17
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% ICaL')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% ICaL')
elseif current_index_SA == 19
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% NCX')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% NCX')
elseif current_index_SA == 21
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% SERCA')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% SERCA')
elseif current_index_SA == 22
    subplot(2,length(output_index_array),1), title('Mouse-to-Human, 50% RyR')
    subplot(2,length(output_index_array),5), title('Rabbit-to-Human, 50% RyR')
end

%% Show summary of translation results
if current_index_SA == 0 % Drug-free
    outputs_mouse = outputs_mouse_baseline;
    outputs_rabbit = outputs_rabbit_baseline;
    outputs_human = outputs_human_baseline;
else
    outputs_mouse = outputs_mouse_block(current_index,:);
    outputs_rabbit = outputs_rabbit_block(current_index,:);
    outputs_human = outputs_human_block(current_index,:);
end

% Simulated data
APD90_mouse = outputs_mouse(5);
APD50_mouse = outputs_mouse(7);
CaTamp_mouse = outputs_mouse(11);
CaTt50_mouse = outputs_mouse(13);
CaTtau_mouse = outputs_mouse(14);

APD90_rabbit = outputs_rabbit(5);
APD50_rabbit = outputs_rabbit(7);
CaTamp_rabbit = outputs_rabbit(11);
CaTt50_rabbit = outputs_rabbit(13);
CaTtau_rabbit = outputs_rabbit(14);

APD90_human = outputs_human(5);
APD50_human = outputs_human(7);
CaTamp_human = outputs_human(11);
CaTt50_human = outputs_human(13);
CaTtau_human = outputs_human(14);

%all_predictions_from_mouse
APD90_from_mouse = [all_predictions_from_mouse(1,4) all_predictions_from_mouse(2,1) all_predictions_from_mouse(3,1) all_predictions_from_mouse(4,1)];
APD50_from_mouse = [all_predictions_from_mouse(1,5) all_predictions_from_mouse(2,2) all_predictions_from_mouse(3,2) all_predictions_from_mouse(4,2)];
CaTt50_from_mouse = [all_predictions_from_mouse(1,9) all_predictions_from_mouse(2,5) all_predictions_from_mouse(3,3) 0];
CaTtau_from_mouse = [all_predictions_from_mouse(1,10) all_predictions_from_mouse(2,6) all_predictions_from_mouse(3,4) 0];
CaTamp_from_mouse = [all_predictions_from_mouse(1,7) all_predictions_from_mouse(2,4) 0 0];

%all_predictions_from_rabbit
APD90_from_rabbit = [all_predictions_from_rabbit(1,4) all_predictions_from_rabbit(2,1) all_predictions_from_rabbit(3,1) all_predictions_from_rabbit(4,1)];
APD50_from_rabbit = [all_predictions_from_rabbit(1,5) all_predictions_from_rabbit(2,2) all_predictions_from_rabbit(3,2) all_predictions_from_rabbit(4,2)];
CaTt50_from_rabbit = [all_predictions_from_rabbit(1,9) all_predictions_from_rabbit(2,5) all_predictions_from_rabbit(3,3) 0];
CaTtau_from_rabbit = [all_predictions_from_rabbit(1,10) all_predictions_from_rabbit(2,6) all_predictions_from_rabbit(3,4) 0];
CaTamp_from_rabbit = [all_predictions_from_rabbit(1,7) all_predictions_from_rabbit(2,4) 0 0];

figure,set(gcf,'color','w','Position',[50,100,1500,750])
subplot(2,5,1),hold on
b = bar([APD90_mouse flip(APD90_from_mouse) APD90_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_mouse; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD90 (ms)')
xticks(1:6); xticklabels({'M','P2','P4','P6','P10','H'})

subplot(2,5,2),hold on
b = bar([APD50_mouse flip(APD50_from_mouse) APD50_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_mouse; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD50 (ms)')
xticks(1:6); xticklabels({'M','P2','P4','P6','P10','H'})

subplot(2,5,3),hold on
b = bar([CaTt50_mouse flip(CaTt50_from_mouse) CaTt50_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_mouse; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTt50 (ms)')
xticks(1:6); xticklabels({'M','P2','P4','P6','P10','H'})
%title('Mouse-to-Human')

subplot(2,5,4),hold on
b = bar([CaTtau_mouse flip(CaTtau_from_mouse) CaTtau_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_mouse; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTtau (ms)')
xticks(1:6); xticklabels({'M','P2','P4','P6','P10','H'})

subplot(2,5,5),hold on
b = bar(1e6*[CaTamp_mouse flip(CaTamp_from_mouse) CaTamp_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_mouse; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTamp (nM)')
xticks(1:6); xticklabels({'M','P2','P4','P6','P10','H'})

subplot(2,5,6),hold on
b = bar([APD90_rabbit flip(APD90_from_rabbit) APD90_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_rabbit; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD90 (ms)')
xticks(1:6); xticklabels({'R','P2','P4','P6','P10','H'})

subplot(2,5,7),hold on
b = bar([APD50_rabbit flip(APD50_from_rabbit) APD50_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_rabbit; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('APD50 (ms)')
xticks(1:6); xticklabels({'R','P2','P4','P6','P10','H'})

subplot(2,5,8),hold on
b = bar([CaTt50_rabbit flip(CaTt50_from_rabbit) CaTt50_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_rabbit; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTt50 (ms)')
xticks(1:6); xticklabels({'R','P2','P4','P6','P10','H'})
%title('Rabbit-to-Human')

subplot(2,5,9),hold on
b = bar([CaTtau_rabbit flip(CaTtau_from_rabbit) CaTtau_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_rabbit; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTtau (ms)')
xticks(1:6); xticklabels({'R','P2','P4','P6','P10','H'})

subplot(2,5,10),hold on
b = bar(1e6*[CaTamp_rabbit flip(CaTamp_from_rabbit) CaTamp_human],'FaceColor',[1 1 1]);
b.FaceColor = 'flat'; b.CData(1,:) = color_rabbit; b.CData(6,:) = color_human;
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('CaTamp (nM)')
xticks(1:6); xticklabels({'R','P2','P4','P6','P10','H'})

if current_index_SA == 0
    subplot(2,5,3), title('Mouse-to-Human, control')
    subplot(2,5,8), title('Rabbit-to-Human, control')
elseif current_index_SA == 2
    subplot(2,5,3), title('Mouse-to-Human, 50% INaL')
    subplot(2,5,8), title('Rabbit-to-Human, 50% INaL')
elseif current_index_SA == 5
    subplot(2,5,3), title('Mouse-to-Human, 50% Itof')
    subplot(2,5,8), title('Rabbit-to-Human, 50% Itof')
elseif current_index_SA == 8
    subplot(2,5,3), title('Mouse-to-Human, 50% IKr')
    subplot(2,5,8), title('Rabbit-to-Human, 50% IKr')
elseif current_index_SA == 13
    subplot(2,5,3), title('Mouse-to-Human, 50% IK1')
    subplot(2,5,8), title('Rabbit-to-Human, 50% IK1')
elseif current_index_SA == 17
    subplot(2,5,3), title('Mouse-to-Human, 50% ICaL')
    subplot(2,5,8), title('Rabbit-to-Human, 50% ICaL')
elseif current_index_SA == 19
    subplot(2,5,3), title('Mouse-to-Human, 50% NCX')
    subplot(2,5,8), title('Rabbit-to-Human, 50% NCX')
elseif current_index_SA == 21
    subplot(2,5,3), title('Mouse-to-Human, 50% SERCA')
    subplot(2,5,8), title('Rabbit-to-Human, 50% SERCA')
elseif current_index_SA == 22
    subplot(2,5,3), title('Mouse-to-Human, 50% RyR')
    subplot(2,5,8), title('Rabbit-to-Human, 50% RyR')
end