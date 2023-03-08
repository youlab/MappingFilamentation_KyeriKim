clear; close all;
%% Basal set up
% shuffle the random number generator
rng('shuffle'); % For "rand" later
%% -------------------- Cell length stochastic simulation --------------------

%% Initial set up

% simulation time window
t_window        = 10;  %in hour
time_step       = 0.1; %in min
delta_tau       = time_step/60;             % time interval of time step (tstep) # Kyeri modified
minDataLen      = t_window/delta_tau;       % biomass data collected in this window, #Kyeri reduces it down 
time_vector     = delta_tau*(0:minDataLen);	% generate time steps
t_limit         = time_vector(end);     	% maximum time for data collection\

% cell numbers and sizes
nMax            = 10000;      % maximum number of cells 10000
nI              = nMax;       % if considering division, need this parameter! initial number of cells (nI < nMax)
LI_0            = 4;          % initial cell size

% initial biomass, growth rate matrices: for noise control
biomass=zeros(minDataLen+1,nMax); %row: time steps including t0, column: cell ID 
cell_GR = zeros(1,nMax);

%% PL assign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lysis profile generation:
Lc              = 24;
H               = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cell growth and Noise Level 

% average growth rate
avg_mu          = 0.4;

% Initial cell size noise level
sigma3          = 0.1*LI_0; 

% growth rate noise level
sigma1          = 0.1*avg_mu; 

%% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial cell biomass + noise
biomass(1,:)= LI_0+sigma3*randn(1,nI);

% all cell length should be > 0
for cID=1:nI
    while biomass(1,cID)<=0; 
        biomass(1,cID)=LI_0+sigma3*randn;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial growth rates, + noise
cell_GR = avg_mu+sigma1*randn(1,nI);

% all growth rate > 0
for cID = 1:nI
    while cell_GR(cID)<=0; 
        cell_GR(cID) = avg_mu+sigma1*randn;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cell length simulation : Kyeri's version
% Definition of P_lysis = L^H/(L^H + Lc^H); applied PL in each time step is
% No division mode: only counts p_lysis by length, t-step wise caldulation
% Hazard function based Cell lysis in each step


for cID = 1:nI
    biomass(:,cID)= biomass(1,cID)*exp(cell_GR(cID)*time_vector);
    % lysis probability
    % PL = biomass(:,cID).^H./(biomass(:,cID).^H+Lc^H);

    % P_survival = 1-PL
    P_survival = 1-[biomass(:,cID).^H./(biomass(:,cID).^H+Lc^H)];
    % rho_L = dPL/dL
    rho_lysis = H*Lc^H*biomass(:,cID).^(H-1)./(biomass(:,cID).^H+Lc^H).^2;
    % P_Haz = rho_L/1-PL
    P_Haz = rho_lysis./P_survival;
    %P_Haz*dL (integration during time interval)
    P_Haz_integ = P_Haz.*[diff(biomass(:,cID));0];
    
    for tstep = 2:minDataLen+1 
        if rand < P_Haz_integ(tstep) % draw one random number and compare to P_Haz_integ, to find when it lysed.
            biomass(tstep+1:end,cID) = 0;
            break
        end
    end
end


%% Plots -----------------------------------------------------------------------
%% PL from final distribution
co= ["#629ABE";"#386A88"];

% Collect final lengths
finalL = zeros(1, nI);
for cID = 1:nI
    finalL(cID) = max(biomass(:,cID));
end
Edge = linspace(0,max(finalL),1000);

%% Final lenths and PL
figure(1);clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1,1]);

% Final Lengths
subplot(141)
    h = histogram(finalL, Edge, 'EdgeColor', co(1));
        xlabel('Final Cell Length', 'fontsize', 15)
        ylabel('Count', 'fontsize', 15)
        axis square
        a = gca;  
        a.FontSize = 20;

    
    Lmid = [Edge(1:end-1)+Edge(2:end)]/2;
    PL_final = cumsum(h.Values)/nI
    
% Assigned lysis profile (P_L)    
subplot(142) 
    P_lysis = Lmid.^H./(Lmid.^H + Lc^H);
    plot(Lmid, P_lysis, 'k-'); hold on
    scatter(Lc, 0.5);
        text(Lc+3,0.5, {['Lc = ',num2str(Lc)];  ['H = ',num2str(H)]}, 'fontsize',16)
        xlabel('Cell Length')
        ylabel('P_L Assigned')
        ylim([0 1])
        xlim([0 50])
        xticks([0 25 50])
        axis('square')
        a = gca;  
        a.FontSize = 20;
    
%PL from final lengths
subplot(143) 
    plot(Lmid, PL_final, 'LineWidth',3, 'color', co(1)); hold on %, 'm-')
    % fit to Hill Equation
    [fitresult, gof] = HillFit_K( Lmid, PL_final, LI_0 );
    scatter(fitresult.c, 0.5);
        text(Lc+3,0.5, {[sprintf('Lc = %.2f', fitresult.c)];  [sprintf('H = %.1f', fitresult.n)]}, 'fontsize',16)
        axis('square')
        xlabel('Cell Length', 'fontsize', 15)
        xlim([0 50])
        xticks([0 25 50])
        ylabel('P_L Simulated', 'fontsize', 15)
        yticks([0 0.5 1])
        a = gca;  
        a.FontSize = 20;
        
% Compare PL assigned - PL simulated       
subplot(144)
    plot(P_lysis, PL_final,'LineWidth',1.5, 'color', co(1)); hold on
    plot(P_lysis,P_lysis, 'k--')
       axis('square')
       xlabel('P_L Assigned', 'fontsize', 15)
       ylabel('P_L Simulated', 'fontsize', 15)
       xticks([0 0.5 1]);
       yticks([0 0.5 1])
       a = gca;  
       a.FontSize = 20;
%% Total biomass
figure(2);clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.71,0.6]);

subplot(121) % number of live cells
    plot(time_vector, sum(biomass~=0,2), 'LineWidth',3, 'color', co(1)) 
        axis('square')
        xlabel('Time (h)', 'fontsize', 15)
        ylabel('Cell number')
        xlim([0 time_vector(end)])
        ylim([0 nI])
        yticks([0 0.5 1]*nI)
        a = gca;  
        a.FontSize = 20;    
        
subplot(122) %normalized biomass = Total biomass/N_I
Le = sum(biomass,2)/nI;
   plot(time_vector, Le, 'LineWidth',3, 'color', co(1)) ; hold on
   yline(Le(1), 'k--')
        axis('square')
        xlabel('Time (h)', 'fontsize', 15)
        xlim([0 time_vector(end)])
        ylabel('Total biomass/N_I', 'fontsize', 15)
        a = gca;  
        a.FontSize = 20;
%% Survivor Length distributions 
figure(3);clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1,1]);

for i = 0:3
subplot(1,4,i+1)
    idx = find(time_vector == 2*i);
    surv = biomass(idx,:);
    surv(surv==0) = nan;
    h = histogram(surv, Edge, 'EdgeColor', co(1));

        if i==0
        xlabel('Cell Length', 'fontsize', 15)
        ylabel('Count', 'fontsize', 15)
        end
        xlim([2 70])
        title([num2str(time_vector(idx)), 'h'])
        axis square
        a = gca;  
        a.FontSize = 20;
end

%% Hill function fitting
function [fitresult, gof] = HillFit_K( Lt, PLs, LI_0 )
idx = length(Lt);
        [xData, yData] = prepareCurveData( [0,Lt(1:idx)], [0,PLs(1:idx)] );
        
        ft = fittype( 'x^n/(c^n+x^n)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.DiffMinChange = 0.000001;
        opts.Display = 'Off';
        opts.Robust = 'Bisquare';%'LAR';
        opts.Lower = [LI_0 0.1];
        opts.StartPoint = [LI_0 0.488158820963235];
        opts.Upper = [10*LI_0 20];

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
end