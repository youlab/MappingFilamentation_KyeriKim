clear all; close all;
%% color set up
red = [0.6350    0.0780    0.1840];
org = [0.8500    0.3250    0.0980];
ylw = [0.9290    0.6940    0.1250];
grn = [0.4660    0.6740    0.1880];
sky = [0.3010    0.7450    0.9330];
blu = [0         0.4470    0.7410];
ppl = [0.4940    0.1840    0.5560];
co=[red;ylw;grn;blu;ppl];
set(groot, 'defaultAxesColorOrder',co);
%% Basal set up
% shuffle the random number generator
rng('shuffle'); % For "rand" later

% Save plot figures
Fdate = datestr(now, 'yyyy_mm_dd_HHMM');
%% -------------------- Cell length stochastic simulation --------------------

%% Initial set up

% simulation time window
t_window        = 8;                            %in hour
time_step       = 0.1;                          %in min
delta_tau       = time_step/60;                 % time interval of time step (tstep) # Kyeri modified
minDataLen      = t_window/delta_tau;           % biomass data collected in this window, #Kyeri reduces it down 
time_vector     = delta_tau*(0:minDataLen+1);	% generate time steps
t_limit         = time_vector(end);             % maximum time for data collection\

% cell numbers and sizes
nMax            = 40000;                        % maximum number of cells 40000
nI              = 2000;                         % if considering division, need this parameter! initial number of cells (nI < nMax)
LI_0            = 4;                            % initial cell size

%% PL assign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antibiotic doses for simulations
As  = [0 5 10 20];
% Lysis profile generation:
Lcs = round(398./As+3); 
% using rounded regression parameters from Fig. 2
%   p1=398.3698    p2=2.9373
H   = 5; 

% Initialize the biomass storing array for boxplots
boxbio = nan(nMax,length(As),3);

% Division length Ld: Ld = alpha*Li+beta;
alpha = 0.871;
beta = 2.7;
%% Cell growth and Noise Level 
% average growth rate
avg_mu          = 0.4;

% Initial cell size noise level
sigma3 = 0.1*LI_0;      %-- 10% of initial size *randn

% growth rate noise level
sigma1 = 0.05*avg_mu;    %-- 5% of avg_mu * randn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial cell size generation: Kyeri modified
%---------------------------------------------------------------------------
        biomass_ini(1,1:nI)=LI_0+sigma3*randn(1,nI);

        % all cell length should be > 0
        for cID = 1:nI
            while biomass_ini(1,cID)<=0
                biomass_ini(1,cID)=LI_0+sigma3*randn;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single cell growth rate vector 
%---------------------------------------------------------------------------
    % initial growth rates, + noise

    cell_GR_ini(1:nI) = avg_mu+sigma1*randn(1,nI);
    for cID = 1:nI
        % all growth rate > 0
        while cell_GR_ini(cID)<=0
            cell_GR_ini(cID) = avg_mu+sigma1*randn;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Biomass simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for aa = 1:length(As)
    A = As(aa);
    Lc = Lcs(aa);
    
    % initialization of biomass, growth rate matrices 
    biomass=NaN(length(time_vector),nMax); %row: time steps including t0, column: cell ID 
    cell_GR = zeros(1,nMax);

    % assign individual initial length and growth rate
    biomass(1, 1:nI) = biomass_ini;
    cell_GR(1:nI) = cell_GR_ini;
    
    % Initialize nMax tracking parameter 
    newborn = nI;
%% ------------------------------------------------------------------------------------------
%% Cell length simulation : Kyeri's version

    % First time step: only elongation
    biomass(2,1:nI) = biomass(1,1:nI).*exp(cell_GR(1:nI)*delta_tau);

    % for each cells, walk through the time steps to decide elongation and
    % lysis (or division in antibiotic=0)
    for cID = 1:nMax 
        % for each cell in biomass tracking matrix, find the initial cell
        % length. (need this step due to division, initial length of
        % daughter cells will be found here)
        LIcID = biomass( find(biomass(:,cID)>0,1)  ,cID ); % set the initial length for division decision.
        % division length
        Ld = alpha*LIcID+beta;
        
            % iteration of "second time step to final time steps"
            for ts = 2:size(biomass,1)-1
                if biomass(ts,cID)>0 % calculate only survivors
                    
                    % lysis probability
                    % PL = biomass(ts,cID)^H/(biomass(ts,cID)^H+Lc^H);
                    
                    % P_survival = 1-PL
                    P_survival = 1-(  biomass(ts,cID)^H/(biomass(ts,cID)^H+Lc^H)  );
                    % rho_L = dPL/dL
                    rho_lysis = (H*Lc^H*biomass(ts,cID)^(H-1))/((biomass(ts,cID)^H+Lc^H)^2);
                    % P_Haz = rho_L/1-PL
                    P_Haz = rho_lysis/P_survival;
                    %P_Haz*dL (integration during time interval)
                    P_Haz_integ = P_Haz*(biomass(ts,cID)-biomass(ts-1,cID)); 
                    

                    % lysis?
                    if A~=0 && rand < P_Haz_integ % with antibiotics, draw one random number and compare to P_Haz_integ, to find when it lysed.
                        biomass(ts+1:end,cID) = 0;

                    % division?
                    elseif A==0 && biomass(ts,cID) > Ld % without antibiotics, when L > Ld, force to divide into two equal length cells
                       
                        % mother cell: track in the same cell ID. 
                        biomass(ts+1,cID) = (biomass(ts,cID)/2)*exp(cell_GR(cID)*delta_tau); % division and elongation
                        LIcID = biomass(ts,cID)/2; % assign new initial length for the mother cell
                        
                        % daughter cell: track in the empty column as a new cell, within nMax cells (no more tracking
                        % after nMax of cells in the lifetime)
                        if newborn < nMax 
                            newborn = newborn+1;
                            cell_GR(newborn) = cell_GR(cID); % assign the same growth rate to daughter cell
                            biomass(ts+1,newborn) = biomass(ts+1,cID);
                        end
                        
                    % elongation! (if no lysis or division)
                    else  biomass(ts+1,cID)= biomass(ts,cID)*exp(cell_GR(cID)*delta_tau);
                    end

                end
            end
    end

%% Plots -----------------------------------------------------------------------
% L effective: normalized length
Le = nansum(biomass,2)'/nI;

% Total biomass of survivors 
figure(1);
    plot(time_vector, Le); hold on
        axis('square')
        xlabel('Time (h)')
        xlim([0 time_vector(end)])
        ylabel('Total biomass/N_I')

% For each antibiotic dose, plot total biomass, tracked cell lengths, 2,4,6h survivor lengths         
figure(aa+10);clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1,1]);
        
% total biomass in normalized plot
subplot(231)
    plot(time_vector, Le)
        axis('square')
        xlabel('Time (h)')
        xlim([0 time_vector(end)])
        ylabel('Total biomass/N_I')

% cell track        
subplot(232)
    plot(time_vector, biomass)

    
samp_int = 2; % 2h interval cell length
Edge = linspace(0,50,250);

% 2h survivor    
subplot(234) 
time_idx = (samp_int/delta_tau)*1+1;
h = histogram(biomass(time_idx,(biomass(time_idx,:)>0)), Edge);
h.FaceColor = co(1,:,:);
boxbio(:,aa,1) = biomass(time_idx,:);
title("2h")
% 4h survivor    
subplot(235) 
time_idx = (samp_int/delta_tau)*2+1;
h = histogram(biomass(time_idx,(biomass(time_idx,:)>0)), Edge);
h.FaceColor = co(2,:,:);
boxbio(:,aa,2) = biomass(time_idx,:);
title("4h")
% 4h survivor    
subplot(236) 
time_idx = (samp_int/delta_tau)*3+1;
h = histogram(biomass(time_idx,(biomass(time_idx,:)>0)), Edge);
h.FaceColor = co(3,:,:);
boxbio(:,aa,3) = biomass(time_idx,:);
title("6h")
end

%% Box plot
figure(200);clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.1, 0.4,0.9]);
colors = ["#C0D5E5"; "#93BBD2";"#629ABE";"#386A88"];
boxbio2 = boxbio;
boxbio2(boxbio ==0) = NaN;


subplot(2,1,1)
for i = 1:4
boxchart(i*ones(size(boxbio2(:,i,1))), boxbio2(:,i,1),'BoxFaceColor', colors(i),'MarkerColor', [0 0 0] , 'MarkerStyle','.'); hold on
end

subplot(2,1,2)
for i = 1:4
boxchart(i*ones(size(boxbio2(:,i,3))), boxbio2(:,i,3),'BoxFaceColor', colors(i),'MarkerColor', [0 0 0] ,  'MarkerStyle','.'); hold on
end

for jj = 1:2
    subplot(2,1,jj)
        xlabel("Carb (µg/ml)")
        xticks([1:4])
        xticklabels(As)
        ylabel("Cell length (µm)")
        xlim([0.5 4.5])
        ylim([0 max(boxbio, [], 'all')+5])
        aaa = gca;
        aaa.FontSize =20
        axis square
end
