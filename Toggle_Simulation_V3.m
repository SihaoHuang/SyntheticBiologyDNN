%% Toggle Setup

%% Biosystem setup
neuron = BioSystem();

% Degradation rates
neuron.AddConstant(Const('k_degrep1', 0.15));
neuron.AddConstant(Const('k_degrep2', 0.15));
neuron.AddConstant(Const('k_degActC', 0.15));

% Hill parameters
neuron.AddConstant(Const('K_rep1', 9));
neuron.AddConstant(Const('K_rep2', 9));
neuron.AddConstant(Const('K_IndT', 1));
neuron.AddConstant(Const('K_ActC', 1));
neuron.AddConstant(Const('n_rep1', 9));
neuron.AddConstant(Const('n_rep2', 9));
neuron.AddConstant(Const('n_IndT', 1));
neuron.AddConstant(Const('n_ActC', 1));

% Production rates
neuron.AddConstant(Const('k_seqrep1', 0.02));
neuron.AddConstant(Const('k_seqrep2', 0.02));
neuron.AddConstant(Const('k_prodE', 2));
neuron.AddConstant(Const('k_prodC', 2));
neuron.AddConstant(Const('k_prodD', 0.2));

% Define compositors
drep1dt = neuron.AddCompositor("rep1", 0);
drep2dt = neuron.AddCompositor("rep2", 0);
dInd1dt = neuron.AddCompositor("Ind1", 0.1);
dInd2dt = neuron.AddCompositor("Ind2", 0.1);
dIndTdt = neuron.AddCompositor("IndT", 0);
dActCdt = neuron.AddCompositor("ActC", 0);

% Inputs (assumed to be constant and controlled)
neuron.AddPart(Part('Ind1', [dInd1dt], [Rate('0')]));
neuron.AddPart(Part('Ind2', [dInd2dt], [Rate('0')]));
neuron.AddPart(Part('IndT', [dIndTdt], [Rate('0')]));

% Dynamical equations
% Strand C, both repressions are OR'd
neuron.AddPart(Part('rep1', [drep1dt], [Rate('k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) - k_seqrep1 * rep1 * Ind1 - k_degrep1 * rep1')]));

% Strand E, ctivation and repression are OR'd
neuron.AddPart(Part('rep2', [drep2dt], [Rate('k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) + k_prodC * (ActC ^ n_ActC / ((K_ActC ^ n_ActC) + ActC ^ n_ActC)) - k_seqrep2 * rep2 * Ind2 - k_degrep2 * rep2')]));

% Strand D
neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));
% neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));
%% Time domain 

[T, Y] = neuron.run([0 100])

figure();

subplot(2, 1, 1);
plot(T, Y(:, neuron.CompositorIndex('rep1')))
ylim([0 20]);
legend('rep1')
xlabel('Time (s)');
ylabel('Concentration (nM)');

subplot(2, 1, 2);
plot(T, Y(:, neuron.CompositorIndex('rep2')))
ylim([0 20]);
legend('rep2')
xlabel('Time (s)');
ylabel('Concentration (nM)');

%% Sweep
% 
% figure()
% ax = gca();
% hold(ax, 'on')
% ax.ColorOrder = jet(20);
% 
% s = 100;
% out = zeros(s);
% x_axis = zeros(s);
% 
% for t = 1:16
%     neuron.ChangeInitialValue("IndT", t/500);
%     for x = 1:s
%         neuron.ChangeInitialValue("Ind1", x/100);
%         [T, Y] = neuron.run([0 1000]);
%         out(x) = Y(25, neuron.CompositorIndex('rep2'));
%         x_axis(x) = x/100;
%     end
%     ax.ColorOrderIndex = t
%     fprintf('index is %d\n',ax.ColorOrderIndex);
%     plot(ax, x_axis, out);
% end
% 
% legend('[IndT] = 0.002 nM', '[IndT] = 0.004nM', '[IndT] = 0.006nM', '[IndT] = 0.008nM', '[IndT] = 0.010nM', '[IndT] = 0.012nM', '[IndT] = 0.014nM', '[IndT] = 0.016nM', '[IndT] = 0.018nM', '[IndT] = 0.020nM', ...
%     '[IndT] = 0.022nM', '[IndT] = 0.024nM', '[IndT] = 0.026nM', '[IndT] = 0.028nM', '[IndT] = 0.030nM', '[IndT] = 0.032nM');
% 
% xlabel('Ind1');
% ylabel('rep1 Concentration (nM)');