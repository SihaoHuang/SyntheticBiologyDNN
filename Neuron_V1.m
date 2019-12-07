% Neuron configuration
input1 = 0;
input2 = 0;
weight1 = 0; 
weight2 = 0;
threshold = 0; 

%% Biosystem setup
neuron = BioSystem();

% Degradation rates
neuron.AddConstant(Const('k_degrep1', 0.15));
neuron.AddConstant(Const('k_degrep2', 0.15));
neuron.AddConstant(Const('k_degActC', 0.5));
neuron.AddConstant(Const('k_degActE', 0.5));
neuron.AddConstant(Const('k_degir1', 0.5));
neuron.AddConstant(Const('k_degout', 0.5));
neuron.AddConstant(Const('k_degLacI', 0.5));
neuron.AddConstant(Const('k_degTetR', 0.5));
% Input-related parameters
neuron.AddConstant(Const('K_A', 5));
neuron.AddConstant(Const('K_B', 5));
neuron.AddConstant(Const('n_A', 1));
neuron.AddConstant(Const('n_B', 1));
neuron.AddConstant(Const('k_deg', 2));
% Hill parameters
neuron.AddConstant(Const('K_LacI', 2));
neuron.AddConstant(Const('K_TetR', 2));
neuron.AddConstant(Const('K_IndT', 2));
neuron.AddConstant(Const('K_ActC', 2));
neuron.AddConstant(Const('K_rep1', 9));
neuron.AddConstant(Const('K_ActE', 2));
neuron.AddConstant(Const('K_rep2', 9));
neuron.AddConstant(Const('n_LacI', 2));
neuron.AddConstant(Const('n_TetR', 2));
neuron.AddConstant(Const('n_IndT', 1));
neuron.AddConstant(Const('n_ActC', 2));
neuron.AddConstant(Const('n_rep1', 5));
neuron.AddConstant(Const('n_ActE', 2));
neuron.AddConstant(Const('n_rep2', 5));
% Production rates
neuron.AddConstant(Const('k_seqLacI', 100));
neuron.AddConstant(Const('k_seqTetR', 100));
neuron.AddConstant(Const('k_seqrep1', 100));
neuron.AddConstant(Const('k_prodA', 2));
neuron.AddConstant(Const('k_prodB', 2));
neuron.AddConstant(Const('k_prodC', 2));
neuron.AddConstant(Const('k_prodD', 10));
neuron.AddConstant(Const('k_prodE', 2));
neuron.AddConstant(Const('k_prodLacI', 10));
neuron.AddConstant(Const('k_prodTetR', 10));

% Define compositors
dLacIdt = neuron.AddCompositor("LacI", 0);
dTetRdt = neuron.AddCompositor("TetR", 0);
drep1dt = neuron.AddCompositor("rep1", 0);
dir1dt = neuron.AddCompositor("ir1", 0);
dActCdt = neuron.AddCompositor("ActC", 0);
dActEdt = neuron.AddCompositor("ActE", 0);
drep2dt = neuron.AddCompositor("rep2", 0);
doutdt = neuron.AddCompositor("out", 0);

dIPTGdt = neuron.AddCompositor("IPTG", input1);
daTcdt = neuron.AddCompositor("aTc", input2);
dIndTdt = neuron.AddCompositor("IndT", threshold);
dfAdt = neuron.AddCompositor("fA", weight1);
dfBdt = neuron.AddCompositor("fB", weight2);

% Inputs (assumed to be constant and controlled)
neuron.AddPart(Part('IPTG', [dIPTGdt], [Rate('0')]));
neuron.AddPart(Part('aTc', [daTcdt], [Rate('0')]));
neuron.AddPart(Part('IndT', [dIndTdt], [Rate('0')]));
neuron.AddPart(Part('fA', [dfAdt], [Rate('0')]));
neuron.AddPart(Part('fB', [dfBdt], [Rate('0')]));

% Dynamical equations
neuron.AddPart(Part('LacI', [dLacIdt], [Rate('k_prodLacI - k_degLacI * LacI - (fA/(fA - K_A))*k_deg*LacI - k_seqLacI * LacI * IPTG')]));

neuron.AddPart(Part('TetR', [dTetRdt], [Rate('k_prodTetR - k_degTetR * TetR - (fB/(fB - K_B))*k_deg*TetR - k_seqTetR * TetR * aTc')]));

neuron.AddPart(Part('ActE', [dActEdt], [Rate('k_prodA * (K_LacI ^ n_LacI)/((K_LacI ^ n_LacI) + LacI ^ n_LacI) + k_prodB * (K_TetR ^ n_TetR)/((K_TetR ^ n_TetR) + TetR ^ n_TetR) - k_degActE * ActE')]));

neuron.AddPart(Part('rep1', [drep1dt], [Rate('k_prodA * (K_LacI ^ n_LacI)/((K_LacI ^ n_LacI) + LacI ^ n_LacI) + k_prodB * (K_TetR ^ n_TetR)/((K_TetR ^ n_TetR) + TetR ^ n_TetR) + k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) * (ActE ^ n_ActE / ((K_ActE ^ n_ActE) + ActE ^ n_ActE)) - k_seqrep1 * rep1 * ir1 - k_degrep1 * rep1')]));
    
neuron.AddPart(Part('out', [doutdt], [Rate('k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) * (ActE ^ n_ActE / ((K_ActE ^ n_ActE) + ActE ^ n_ActE)) - k_degout * out')]));

neuron.AddPart(Part('rep2', [drep2dt], [Rate('k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) * (ActC ^ n_ActC / ((K_ActC ^ n_ActC) + ActC ^ n_ActC)) + k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degrep2 * rep2')]));

neuron.AddPart(Part('ir1', [dir1dt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degir1 * ir1')]));

neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));

%% Run the system
[T, Y] = neuron.run([0 100])

% plot
figure();

subplot(4, 1, 1);
plot(T, Y(:, neuron.CompositorIndex('out')))
ylim([0 20]);
legend('out')
xlabel('Time (s)');
ylabel('Concentration (nM)');

subplot(4, 1, 2);
plot(T, Y(:, neuron.CompositorIndex('rep1')), ...
     T, Y(:, neuron.CompositorIndex('rep2')))
% ylim([0 50]);
legend('rep1', 'rep2')
xlabel('Time (s)');
ylabel('Concentration (nM)');

subplot(4, 1, 3);
plot(T, Y(:, neuron.CompositorIndex('ActC')), ...
     T, Y(:, neuron.CompositorIndex('ActE')))
% ylim([0 50]);
legend('ActC', 'ActE')
xlabel('Time (s)');
ylabel('Concentration (nM)');

subplot(4, 1, 4);
plot(T, Y(:, neuron.CompositorIndex('LacI')), ...
     T, Y(:, neuron.CompositorIndex('TetR')))
% ylim([0 12]);
legend('LacI', 'TetR')
xlabel('Time (s)');
ylabel('Concentration (nM)');