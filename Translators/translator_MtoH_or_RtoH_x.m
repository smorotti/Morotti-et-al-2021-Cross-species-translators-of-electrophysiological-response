%% Build Mouse-to-Human or Rabbit-to-Human translator using different sets
% of simulated features

clear
close all
clc

%% Selection of input species
% Set "input_species" to 1 or 2 to select Mouse or Rabbit as input species
input_species = 1; % 1 for Mouse, 2 for Rabbit

%% Selection of considered features
% Change "output_index" from 0 to 4 to select the desired set of features
output_index = 3;

% List of all features available:
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp

if output_index == 0 % ALL
    output_selection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]; % nf = 18
    % ALL OUTPUTS
elseif output_index == 1
    output_selection = [1 3 4 5 7 10 11 12 13 14]; % nf = 10
    % UV	MDP	APamp	APD90	APD50	CaTmin	CaTamp	CaTttp	CaTt50	CaTtau
elseif output_index == 2
	output_selection = [5 7 10 11 13 14]; % nf = 6
    % APD90	APD50	CaTmin	CaTamp	CaTt50	CaTtau
elseif output_index == 3
	output_selection = [5 7 13 14]; % nf = 4
    % APD90	APD50	CaTt50	CaTtau
elseif output_index == 4
	output_selection = [5 7]; % nf = 2
    % APD90	APD50
end

%% Load simulated data
disp('1-Hz');

load outputs_matrix_1500_0p1_300s_mouse_COMMON
all_outputs_mouse_control = all_outputs;

load outputs_matrix_1500_0p1_300s_rabbit_COMMON
all_outputs_rabbit_control = all_outputs;

load outputs_matrix_1500_0p1_300s_human_COMMON
all_outputs_human_control = all_outputs;

if input_species == 1
    disp('Mouse to Human')
    Xblock = all_outputs_mouse_control;
    Yblock = all_outputs_human_control;
elseif input_species == 2   
    disp('Rabbit to Human')
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
 
%% Regression

% Call the PLS routine 
disp('Fitting')
N_fitting = good_count;

% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% N_pars, N_outputs, N_trials, good_count

%% Goodness of fit
% R^2
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

%% Plots
% Plotting options
color = [0 0 0]; % black

plot_matrix = 1;
plot_fitting = 0;
plot_validation = 1;
plot_residuals = 0;

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

if plot_residuals == 1
    % The error terms are assumed to be:
    % 1) Normally distribuited
    % 2) Homoscedatic (same variance at every X)
    % 3) Independent
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
                plot(factor*exp(Yhat(:,dex1)),(exp(Y(:,dex1))-exp(Yhat(:,dex1)))/rSD(dex1),'Marker','o','LineStyle','none','Color',color);
                xlabel(['Predicted ', output_names_Y{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD(dex1)/oSD(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',12)
                % Plot identity line
                xlim_ind = get(gca,'xlim') ;
                plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k') 
                xlim([xlim_ind(1), xlim_ind(2)])
                
                dex1 = dex1+1;
            end
        end
    end
    
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    subplot(2,2,1),bar(R2_fit,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels( gca(), 90)
    title('R^2 values')

    subplot(2,2,2),bar(rSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD./oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')
end

%% Validation
disp('Validation')
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

if plot_residuals == 1
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
                plot(factor*predicted_outputs(:,dex1),(actual_outputs(:,dex1)-predicted_outputs(:,dex1))/rSD_val(dex1),'Marker','o','LineStyle','none');
                xlabel(['Predicted ', output_names_Y{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD_val(dex1)/oSD_val(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',14)
                % Plot identity line
                xlim_ind = get(gca,'xlim') ;
                plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k') 
                xlim([xlim_ind(1), xlim_ind(2)])
                
                dex1 = dex1+1;
            end
        end
    end
    
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    subplot(2,2,1),bar(R2)%,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels(gca(), 90)
    title('R^2 values')

    % Add residuals also in this population
    subplot(2,2,2),bar(rSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD_val./oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')
end
            
n_pop_fitting = N_fitting%;
avg_R2_fitting = avg_R2_fit%;
avg_R2_validation = avg_R2_val%;

%% Set plot_correlation to 1 to investigate correlation among the features
% in input species (1st figure) or output species (2nf figure)
plot_correlation = 0;

% X - input species
index_pci = 1;
if plot_correlation == 1
    figure,set(gcf,'color','w')
    for i = 1:N_outputs_X
        for j = 1:N_outputs_X
            subplot(N_outputs_X,N_outputs_X,index_pci)
            set(gca,'box','off','tickdir','out','fontsize',10)
            index_pci = index_pci+1;
            plot(good_outputs_X(:,j),good_outputs_X(:,i),'k.')
            if i == N_outputs_X
                xlabel(output_names{j})
            end
            if j == 1
                ylabel(output_names{i})
            end
            if i == j
                set(gca,'Color',[0.5 0.5 0.5])
            end
                
            % Correlation
            % Rule of Thumb for Interpreting the Size of a Correlation Coefficient
            % Size of Correlation   Interpretation
            % .90 to 1.00    		Very high positive (negative) correlation
            % .70 to .90     		High positive (negative) correlation
            % .50 to .70     		Moderate positive (negative) correlation
            % .30 to .50     		Low positive (negative) correlation
            % .00 to .30     		Negligible correlation
            [R, P] = corrcoef(good_outputs_X(:,j),good_outputs_X(:,i));
            %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
            title(['R^2 = ',num2str(R(1,2)^2,3)]);
        end
    end
end

% Y - output species
index_pci = 1;
if plot_correlation == 1
    figure,set(gcf,'color','w')
    for i = 1:N_outputs_Y
        for j = 1:N_outputs_Y
            subplot(N_outputs_Y,N_outputs_Y,index_pci)
            set(gca,'box','off','tickdir','out','fontsize',10)
            index_pci = index_pci+1;
            plot(good_outputs_Y(:,j),good_outputs_Y(:,i),'.')
            if i == N_outputs_Y
                xlabel(output_names{j})
            end
            if j == 1
                ylabel(output_names{i})
            end
            if i == j
                set(gca,'Color',[0.5 0.5 0.5])
            end
                
            % Correlation
            % Rule of Thumb for Interpreting the Size of a Correlation Coefficient
            % Size of Correlation   Interpretation
            % .90 to 1.00    		Very high positive (negative) correlation
            % .70 to .90     		High positive (negative) correlation
            % .50 to .70     		Moderate positive (negative) correlation
            % .30 to .50     		Low positive (negative) correlation
            % .00 to .30     		Negligible correlation
            [R, P] = corrcoef(good_outputs_Y(:,j),good_outputs_Y(:,i));
            %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
            title(['R^2 = ',num2str(R(1,2)^2,3)]);
        end
    end
end
