close all
clear all

global Lc H mu tstep
L0  = 4;
mus = 0.5;
Lcs = [4:0.05:20,20.2:0.2:50];
H   = 3;
n0  = 10^7;
%% Colormap setup
cmap_custom = zeros(3, 256);
c0 = 0.9*[1,1,1];
c1 = [0.2081, 0.1663, 0.5292];
c2 = [0.0265, 0.6137, 0.8135];
c3 = [0.1986, 0.7214, 0.6310];
c4 = [0.9856, 0.7372, 0.2537];

for i = 1:39
    cmap_custom(:,i) = c1 + (c0-c1)*i/39;
end
for i = 1:40
    cmap_custom(:,39+i) = c2 + (c0-c2)*(40-i)/(40);
end
for i = 1:40
    cmap_custom(:,79+i) = c2 + (c3-c2)*i/(40);
end
for i = 1:(256-119)
    cmap_custom(:,119+i) = c3 + (c4-c3)*i/(256-119);
end
%% Total biomass simulations

tstep   = 0.5/60/mus; %in min
twindow = 8;
tspan   = 0:tstep:twindow;

biomass = zeros(length(Lcs), length(tspan),length(mus));
%%
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        for i = 1:length(Lcs)
            Lc = Lcs(i);
        for m = 1:length(mus);
            mu = mus(m);    
        [t, y] = ode45(@odelength, tspan, [L0, n0], options);
        % Nt = y(:,2);
        % Lt = y(:,1); 
        biomass(i, :, m) = y(:,2).*y(:,1)/n0;
        end
        end
%% Increased biomass: subtract initial
biomass2 = biomass-biomass(1,1,1); 
%% Plots
figure(1);clf;
    for i = 1:length(Lcs);
        Lc = Lcs(i);
    for m = 1:length(mus);
        mu = mus(m);  
    subplot(1,length(mus),m)
        scatter((mus(m)*tspan), log((Lcs(i)/L0)^H*ones(size(tspan))), 12,biomass2(i,:,m), 'filled'); hold on
    end
    end
%% Analytical Boundary 
tspan2 = 0:0.1:twindow;
phi_all = (   (exp(mus'*tspan2)).^H -exp(mus'*tspan2)  )./(exp(mus'*tspan2)-1);
%%
figure(1)
for m = 1:length(mus);
    mu = mus(m);
subplot(1,length(mus),m)
plot((mus(m)*tspan2), log(phi_all(m,:)),'r-' ,'linewidth',3);
box on
end
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.4,0.7]);
        xlabel('$\mu t$','Interpreter','latex','FontSize',15)
        xlim([0 4])
        ylabel('$\ln{(\frac{L_C}{L_0})^H}$','Interpreter','latex','FontSize',15)
        yticks([1:2:7])
        ylim([1 7])
        axis square
        a = gca;  
        set(gca,'fontsize',30)  
        title('Increased Biomass', 'FontSize',25)
        colormap(cmap_custom')
        colorbar 
%%
figure(1)
axis square
set(gca, 'Position', [-.1 .1 1.2 0.85])
set(gca, 'OuterPosition', [0 0 1 1])
%% ODE
function dydt = odelength(t,y)
global Lc H mu tstep 
L=y(1);
n=y(2);
PL = L^H/(L^H+Lc^H);
pL = H*Lc^H*L^(H-1)/(L^H+Lc^H)^2;
HL = pL/(1-PL)*mu*L;%*tstep;

dLdt = mu*L;
dndt = -HL*n;

dydt=[dLdt;dndt];
end
