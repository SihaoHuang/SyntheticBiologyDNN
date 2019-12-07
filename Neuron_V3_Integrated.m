%% Toggle Setup

%% Biosystem setup
neuron = BioSystem();

% Weights
neuron.AddConstant(Const('K_A', 5));
neuron.AddConstant(Const('K_B', 5));
neuron.AddConstant(Const('n_A', 1));
neuron.AddConstant(Const('n_B', 1));
neuron.AddConstant(Const('k_deg', 2));
neuron.AddConstant(Const('k_seqLacI', 0.1));
neuron.AddConstant(Const('k_seqTetR', 0.1));

% General degradation rates
neuron.AddConstant(Const('k_degrep1', 0.15));
neuron.AddConstant(Const('k_degrep2', 0.15));
neuron.AddConstant(Const('k_degActC', 0.15));
neuron.AddConstant(Const('k_degInd1', 0.15));
neuron.AddConstant(Const('k_degLacI', 0.15));
neuron.AddConstant(Const('k_degTetR', 0.15));
neuron.AddConstant(Const('k_degInd2', 0.15));

% General Hill parameters
neuron.AddConstant(Const('K_rep1', 9));
neuron.AddConstant(Const('K_rep2', 9));
neuron.AddConstant(Const('K_IndT', 1));
neuron.AddConstant(Const('K_ActC', 1));
neuron.AddConstant(Const('n_rep1', 9));
neuron.AddConstant(Const('n_rep2', 9));
neuron.AddConstant(Const('n_IndT', 1));
neuron.AddConstant(Const('n_ActC', 1));
neuron.AddConstant(Const('K_LacI', 2));
neuron.AddConstant(Const('K_TetR', 2));
neuron.AddConstant(Const('n_LacI', 2));
neuron.AddConstant(Const('n_TetR', 2));

% Production rates
neuron.AddConstant(Const('k_prodA', 2));
neuron.AddConstant(Const('k_prodB', 2));
neuron.AddConstant(Const('k_prodC', 2));
neuron.AddConstant(Const('k_prodD', 0.2));
neuron.AddConstant(Const('k_prodE', 2));
neuron.AddConstant(Const('k_prodLacI', 2));
neuron.AddConstant(Const('k_prodTetR', 2));
neuron.AddConstant(Const('k_prodInd2', 0.5));

% Toggle sequestration
neuron.AddConstant(Const('k_seqrep1', 0.02));
neuron.AddConstant(Const('k_seqrep2', 0.02));

% Define compositors
drep1dt = neuron.AddCompositor("rep1", 0);
drep2dt = neuron.AddCompositor("rep2", 0);
dInd1dt = neuron.AddCompositor("Ind1", 0);
dInd2dt = neuron.AddCompositor("Ind2", 0);
dIndTdt = neuron.AddCompositor("IndT", 0);
dActCdt = neuron.AddCompositor("ActC", 0);
dfAdt = neuron.AddCompositor("fA", 0); % weight 1
dfBdt = neuron.AddCompositor("fB", 0); % weight 2
dLacIdt = neuron.AddCompositor("LacI", 0);
dTetRdt = neuron.AddCompositor("TetR", 0);

dIPTGdt = neuron.AddCompositor("IPTG", 0);
daTcdt = neuron.AddCompositor("aTc", 0);

% Inputs (assumed to be constant and controlled)
neuron.AddPart(Part('IndT', [dIndTdt], [Rate('0')]));
neuron.AddPart(Part('fA', [dfAdt], [Rate('0')]));
neuron.AddPart(Part('fB', [dfBdt], [Rate('0')]));
neuron.AddPart(Part('IPTG', [dIPTGdt], [Rate('0')]));
neuron.AddPart(Part('aTc', [daTcdt], [Rate('0')]));

% Dynamical equations

% Weight elements
neuron.AddPart(Part('LacI', [dLacIdt], [Rate('k_prodLacI - k_degLacI * LacI - (fA/(fA - K_A))*k_deg*LacI - k_seqLacI * LacI * IPTG')]));
neuron.AddPart(Part('TetR', [dTetRdt], [Rate('k_prodTetR - k_degTetR * TetR - (fB/(fB - K_B))*k_deg*TetR - k_seqTetR * TetR * aTc')]));

% Connection between inputs and toggle
% Ind1 is under LacI and TetR control
neuron.AddPart(Part('Ind1', [dInd1dt], [Rate('k_prodA * (K_LacI ^ n_LacI)/((K_LacI ^ n_LacI) + LacI ^ n_LacI) + k_prodB * (K_TetR ^ n_TetR)/((K_TetR ^ n_TetR) + TetR ^ n_TetR) - k_degInd1 * Ind1')]));
% Ind2 is constitutively expressed 
neuron.AddPart(Part('Ind2', [dInd2dt], [Rate('k_prodInd2 - k_degInd2 * Ind2')]));

% Strand C, both repressions are OR'd
neuron.AddPart(Part('rep1', [drep1dt], [Rate('k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) - k_seqrep1 * rep1 * Ind1 - k_degrep1 * rep1')]));

% Strand E, ctivation and repression are OR'd
neuron.AddPart(Part('rep2', [drep2dt], [Rate('k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) + k_prodC * (ActC ^ n_ActC / ((K_ActC ^ n_ActC) + ActC ^ n_ActC)) - k_seqrep2 * rep2 * Ind2 - k_degrep2 * rep2')]));

% Strand D
neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));
% neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));

%% Time domain 

% neuron.ChangeInitialValue("IPTG", 0);
% [T, Y] = neuron.run([0 500]);
% 
% figure();
% 
% subplot(4, 1, 1);
% plot(T, Y(:, neuron.CompositorIndex('rep1')))
% ylim([0 20]);
% legend('rep1')
% xlabel('Time (s)');
% ylabel('Concentration (nM)');
% 
% subplot(4, 1, 2);
% plot(T, Y(:, neuron.CompositorIndex('rep2')))
% ylim([0 20]);
% legend('rep2')
% xlabel('Time (s)');
% ylabel('Concentration (nM)');
% 
% subplot(4, 1, 3);
% plot(T, Y(:, neuron.CompositorIndex('LacI')))
% ylim([0 20]);
% legend('LacI')
% xlabel('Time (s)');
% ylabel('Concentration (nM)');
% 
% subplot(4, 1, 4);
% plot(T, Y(:, neuron.CompositorIndex('Ind1')), ...
%      T, Y(:, neuron.CompositorIndex('Ind2')))
% ylim([0 20]);
% legend('Ind1', 'Ind2');
% xlabel('Time (s)');
% ylabel('Concentration (nM)');

% Currently, Ind1 puts the neuron in the high state (note rep2 is the
% output) due to a hazard at the beginning. The neuron has no means of
% resetting back into the low state. Maybe add Ind2 expression to LacI and
% TetR parts?

%% Transfer function

% figure()
% ax = gca();
% hold(ax, 'on')
% ax.ColorOrder = jet(20);
% 
% s = 100;
% out = zeros(s);
% x_axis = zeros(s);
% 
% for x = 1:s
%     neuron.ChangeInitialValue("IPTG", x/100);
%     [T, Y] = neuron.run([0 1000]);
%     out(x) = Y(length(T), neuron.CompositorIndex('rep2'));
%     x_axis(x) = x/100;
% end
% 
% plot(ax, x_axis, out);


%% Threshold sweep

% figure()
% ax = gca();
% hold(ax, 'on')
% ax.ColorOrder = jet(20);
% 
% s = 100;
% out = zeros(s);
% x_axis = zeros(s);
% 
% for t = 1:20
%     neuron.ChangeInitialValue("IndT", t/10000);
%     for x = 1:s
%         neuron.ChangeInitialValue("IPTG", x/100);
%         [T, Y] = neuron.run([0 1000]);
%         out(x) = Y(length(T), neuron.CompositorIndex('rep2'));
%         x_axis(x) = x/99;
%     end
%     ax.ColorOrderIndex = t
%     fprintf('index is %d\n',ax.ColorOrderIndex);
%     plot(ax, x_axis, out);
% end
% 
% xlabel('Ind1');
% ylabel('rep1 Concentration (nM)');

%% 3D Plot of the Neuron Function

s = 10;
out = zeros(s);

for x = 1:s
    for y = 1:s
        neuron.ChangeInitialValue("IPTG", x);
        neuron.ChangeInitialValue("aTc", y);
        [T, Y] = neuron.run([0 100])
        out(x,y) = Y(50, neuron.CompositorIndex('out'))
    end
end

% surface in 3D
figure;
surf(out,'EdgeColor','None');

