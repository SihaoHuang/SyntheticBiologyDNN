%% Toggle Setup

%% Biosystem setup
neuron = BioSystem();

% Degradation rates
neuron.AddConstant(Const('k_degrep1', 0.04));
neuron.AddConstant(Const('k_degrep2', 0.04));
% Hill parameters
neuron.AddConstant(Const('K_rep1', 7));
neuron.AddConstant(Const('K_rep2', 7));
neuron.AddConstant(Const('n_rep1', 4));
neuron.AddConstant(Const('n_rep2', 4));
% Production rates
neuron.AddConstant(Const('k_seqrep1', 0.01));
neuron.AddConstant(Const('k_seqrep2', 0.01));
neuron.AddConstant(Const('k_prodE', 0.5));
neuron.AddConstant(Const('k_prodC', 0.5));

% Define compositors
drep1dt = neuron.AddCompositor("rep1", 0.01);
drep2dt = neuron.AddCompositor("rep2", 0);
dInd1dt = neuron.AddCompositor("Ind1", 0.03);
dInd2dt = neuron.AddCompositor("Ind2", 0);

% Inputs (assumed to be constant and controlled)
neuron.AddPart(Part('Ind1', [dInd1dt], [Rate('0')]));
neuron.AddPart(Part('Ind2', [dInd2dt], [Rate('0')]));

% Dynamical equations
neuron.AddPart(Part('rep1', [drep1dt], [Rate('k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) - k_seqrep1 * rep1 * Ind1 - k_degrep1 * rep1')]));
neuron.AddPart(Part('rep2', [drep2dt], [Rate('k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) - k_seqrep2 * rep2 * Ind2 - k_degrep2 * rep2')]));

%% Time domain 

[T, Y] = neuron.run([0 500])

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

% Sweep
% s = 100;
% out = zeros(s);
% x_axis = zeros(s);
% 
% for x = 1:s
%     neuron.ChangeInitialValue("Ind1", x/100);
%     [T, Y] = neuron.run([0 1000]);
%     out(x) = Y(40, neuron.CompositorIndex('rep2'));
%     x_axis(x) = x/100;
% end
% 
% plot(x_axis, out);
% xlabel('Ind1');
% ylabel('rep1 Concentration (nM)');