% https://math.stackexchange.com/questions/680852/how-to-plot-a-phase-portrait-for-this-system-of-differential-equations

%% Declare constants
%k_prodE = 2;
%k_prodC = 2;
%K_rep2 = 11;
%K_rep1 = 11;
%n_rep2 = 9;
%n_rep1 = 9;
%k_degrep1 = 0.15;
%k_degrep2 = 0.15;

k_prodE = 0.5;
k_prodC = 0.5;
K_rep2 = 7;
K_rep1 = 7;
n_rep2 = 4;
n_rep1 = 4;
k_degrep1 = 0.04;
k_degrep2 = 0.04;

%% Plot nullclines
% syms x y;
% fimplicit(k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + (y) ^ n_rep2)) - k_degrep1 * (x) == 0);
% hold on
% fimplicit(k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + (x) ^ n_rep1)) - k_degrep2 * (y) == 0);
% xlim([0 20]);
% ylim([0 20]);
% xlabel('rep1');
% ylabel('rep2');
% hold on


%% Plot isoclines
syms x y;
s = 10;
% for c = -50:50
%     c = c/3;
%     if c == 0
%         fimplicit(k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + (y) ^ n_rep2)) - k_degrep1 * (x) == c, 'LineWidth',2);
%         hold on
%         fimplicit(k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + (x) ^ n_rep1)) - k_degrep2 * (y) == c, 'LineWidth',2);
%     end
%     fimplicit(k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + (y) ^ n_rep2)) - k_degrep1 * (x) == c, 'Color',[0,0.7,abs(c/30)]);
%     hold on
%     fimplicit(k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + (x) ^ n_rep1)) - k_degrep2 * (y) == c, 'Color',[1,0.5,abs(c/30)]);
% end

fimplicit(k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + (y) ^ n_rep2)) - k_degrep1 * (x) == 0, 'LineWidth',2);
hold on
fimplicit(k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + (x) ^ n_rep1)) - k_degrep2 * (y) == 0, 'LineWidth',2);

xlim([0 20]);
ylim([0 20]);
xlabel('[rep1]');
ylabel('[rep2]');

%% Plot vector field
% xdom = linspace(0,2,31);
% ydom = linspace(0,2,31);
% 
% [rep1,rep2] = meshgrid(xdom,ydom)
% 
% U = k_prodE * ((K_rep2 ^ n_rep2)./((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) - k_degrep1 * rep1; % drep1/dt
% V = k_prodC * ((K_rep1 ^ n_rep1)./((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) - k_degrep2 * rep2; % drep2/dt
% q = quiverc(rep1,rep2,U,V);
% 
