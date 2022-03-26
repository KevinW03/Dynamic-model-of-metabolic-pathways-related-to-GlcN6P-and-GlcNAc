%% Dynamic changes of all intermediate metabolites simulated by the model in the metabolic pathway (Liu et al., 2016; Saa and Nielsen, 2017)

function BNZR100_proGlcNAc
% Strain BNZR1.06 with the highcopy synthetic pathway of GlcNAc

close all

x = 9; % Nubmer of metabolites
m = 100; % Nubmer of loops
v = 11; % Nubmer of reactions
GlcN6P_yield = zeros(m,1); % To save the yield of GlcN6P during m simulations
GlcNAc_yield = zeros(m,1); % To save the yield of GlcNAc during m simulations
X_hat = zeros(201,x,m); % For the save of metabolites concentrations in m simulation loops

global par % include Km,vman,ki, n

% Start of m sampling loops
for i = 1:m
    % Parameters initialization
    par.Km = 1.5*rand(1,v)+0.5.*ones(1,v);
    par.vmax = ones(1,v);
    par.n = 4*rand(1);
    par.ki =rand(1);
    par.vmax(1) = 0.5;
    par.vmax(5) = 0.7; % HMP 
    par.vmax(6) = 1.8; % EMP
    par.vmax(7) = 0.5; % PSP 
    par.vmax(11) = 2; % The synthetic pathway of GlcNAc 
    par.vmax(8) = 0; % The catabolic pathway of GlcNAc
    par.vmax(9) = 0; % The catabolic pathway of GlcN6P
    par.vmax(4) = 2; % Overexpression of glmS
    
    % initial metabolite concentrations
    x0 = zeros(1,x);
    x0(1,1)=1;
    
    [time,X] = ode23s(@pathw_fun,[0 20],x0); %integration
    
    GlcN6P_yield(i,1)=X(end,4);
    GlcNAc_yield(i,1)=X(end,9);
    % The yield of GlcN6P and GlcNAc during m simulations
    
    % interpolate results to compute mean and std of results
    time_hat = 0:0.1:20;
    
    for k = 1:x
        X_hat(:,k,i) = interp1(time,X(:,k),time_hat);% Interpolation in ench loop
    end
    
end % End of m sampling loops

% plot simulation results
meansX = zeros(x,201);
stdX = zeros(x,201);

for k = 1:x
    dummy(:,:)  = X_hat(:,k,:);
    meansX(k,:) = mean(dummy');
    stdX(k,:)   = std(dummy');
end

labclr = {'k-' 'y-' 'b-' 'm-' 'r-' 'g-' };
j=0;
for k = [1,5,6,7,9,4]
    j=j+1;
    boundedline(time_hat,(meansX(k,:)),(stdX(k,:)),labclr{j},'alpha');
    hold on
end

j=0;
for k = [1,5,6,7,9,4]
    j=j+1;
    plot(time_hat,(meansX(k,:)),labclr{j},'LineWidth',1.5)
    set(gca, 'FontSize', 14);
    axis([0 20 -0.5 3]);
    hold on;
end
plot(time_hat,(meansX(4,:)),labclr{j},'LineWidth',3)
set(gca, 'FontSize', 18,'LineWidth', 1.5);
xlabel('Time');ylabel('Concentration');
axis([0 20 -0.5 3]);

writematrix( [GlcN6P_yield GlcNAc_yield],'yield.xlsx');
% create an Excel file, containing the yield of GlcN6P and GlcNAc during m simulations
end
% end

%%linear pathway ODEs - calculate rates and mass balances
function dx_dt = pathw_fun(~,x)
global par
% influx
v_1= par.vmax(1);
% Reaction v_2
v_2 = par.vmax(2)*x(1)/(x(1)+par.Km(2));
% Reaction v_3
v_3 = par.vmax(3)*x(2)/(x(2)+par.Km(3));
% Reaction v_4
v_4 = par.vmax(4)*x(3)/(x(3)+par.Km(4));
% Reaction v_5
v_5 = par.vmax(5)*x(2)/(x(2)+par.Km(5));
% Reaction v_6
v_6 = par.vmax(6)*x(3)/(x(3)+par.Km(6));
% Reaction v_7
v_7 = par.vmax(7)*x(4)/(x(4)+par.Km(7));
% Reaction v_8
v_8 = par.vmax(8)*x(9)/(x(9)+par.Km(8));
% Reaction v_9
v_9 = par.vmax(9)*x(4)/(x(4)+par.Km(9));
% Reaction v_10
v_10 = par.vmax(10)*x(8)/(x(8)+par.Km(10));
% Reaction v_11
v_11 = par.vmax(11)*x(4)/(x(4)+par.Km(11));

% mass balance
dx_dt(1,1) = v_1-v_2;
dx_dt(2,1) = v_2-v_3-v_5;
dx_dt(3,1) = v_3-v_4-v_6+v_9;
dx_dt(4,1) = v_4-v_7+v_8-v_9+v_10-v_11;
dx_dt(5,1) = v_5;
dx_dt(6,1) = v_6;
dx_dt(7,1) = v_7;
dx_dt(8,1) = -v_10;
dx_dt(9,1) = v_11-v_8;
end
