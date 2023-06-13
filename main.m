% main.m
% computes solutions and plots them against data

clear; %close all;

%% Specify which model. (Optional: also specify which data to use.)
% paramfile = "wright_pars.xlsx";
% datafile  = "WeltData.xlsx";
% modelsolver = @wright_solve;
% modellabel = "Wright et al., 2020";
% numCycles = 2; % how many cycles to simulate

% paramfile = "clark_pars.xlsx";
% % datafile = "McLachlanData.xlsx";
% modelsolver = @clark_solve;
% modellabel = "Harris Clark et al., 2003";
% numCycles = 1; % how many cycles to simulate

% paramfile = "gavina_pars.xlsx";
% datafile  = "WeltData.xlsx";
% modelsolver = @gavina_solve;
% modellabel = "Gavina et al., 2022";
% numCycles = 2; % how many cycles to simulate

paramfile = "margolskee_pars.xlsx";
datafile  = "WeltData.xlsx";
modelsolver = @margolskee_solve;
modellabel = "Markgolskee & Selgrade, 2011";
numCycles = 2; % how many cycles to simulate

%% specify which variables to plot
s.variableNames = {'LH','FSH','E2','P4'};
s.lh  = 1; % these indices must be in the same order as in s.variableNames
s.fsh = 2;
s.e2  = 3;
s.p4  = 4;

%% compute solutions
% if datafile not specificied, compute solutions only
% else compute solutions with LH peaks informed by data
if ~exist('datafile','var')
    [T,sols,vars_i,pars,period] = readrun(paramfile,numCycles,modelsolver);
else
    [T,sols,vars_i,pars,period] = readrun(paramfile,numCycles,modelsolver,datafile);
end

SOLS{s.lh} = sols(vars_i.lh,:);
SOLS{s.fsh} = sols(vars_i.fsh,:);
SOLS{s.e2} = sols(vars_i.e2,:);
SOLS{s.p4} = sols(vars_i.p4,:);

%% make figure

N = length(s.variableNames);
figure('Position', [0 300 ceil(N/2)*500 ceil(N/2)*500]); 
tstart=T(1); tend=T(end);

% labels for subplots
letters = ('a':'z').';
char = num2cell(letters(1:N));
char = char';
char = strcat('(',char,')');

% if no datafile, plot only model solutions
% else, plot solutions with data
if ~exist('datafile','var')
    % make figure
    for i=1:N
        subplot(ceil(N/2),ceil(N/2),i); plot(T,SOLS{i},'DisplayName',cell2mat(s.variableNames(i)),'Color','k','LineWidth',2); hold on;
        legend; xlim([tstart tend]);
        text(0.025,0.95,char{i},'Units','normalized','Interpreter','latex');
    end
else
    rmse = zeros(1,N);

    % read data
    [TIME,DATA,SEM] = readData(s.variableNames,datafile);

    % make figure
    for i=1:N

        subplot(ceil(N/2),ceil(N/2),i); 
        for j=1:numCycles
            errorbar(period*(j-1)+TIME{i},DATA{i},SEM{i},'k.','HandleVisibility','off','LineWidth',1); hold on;
        end

        plot(T,SOLS{i},'DisplayName',cell2mat(s.variableNames(i)),'Color','k','LineWidth',2); hold on;

        xlim([tstart tend]); ylim([0 inf]);
        xlabel('Time (days)','Interpreter','latex');
        ylabel(cell2mat(s.variableNames(i)),'Interpreter','latex');
    
        %compute RMSE
        est = interp1(T,SOLS{i},TIME{i});
        dat = DATA{i};
        rmse(i) = sqrt(sum(abs(est-dat).^2)/length(dat));
        
        title(['$RMSE = $ ',num2str(rmse(i))],'Interpreter','latex');
        text(0.025,0.95,char{i},'Units','normalized','Interpreter','latex');
    end
end

sgtitle(strcat(modellabel, " (period = ", num2str(period), " days)"));
