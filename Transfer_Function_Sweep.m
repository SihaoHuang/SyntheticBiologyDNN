%% Base model
 
% Ratio between ActE and rep1
ratio = 1
% Input
in = 0

% Neuron configuration
threshold = 0; 

% Biosystem setup
neuron = BioSystem();

% Degradation rates
neuron.AddConstant(Const('k_degrep1', 0.75));
neuron.AddConstant(Const('k_degrep2', 0.75));
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
neuron.AddConstant(Const('K_IndT', 2));
neuron.AddConstant(Const('K_ActC', 2));
neuron.AddConstant(Const('K_rep1', 10));
neuron.AddConstant(Const('K_ActE', 2));
neuron.AddConstant(Const('K_rep2', 10));
neuron.AddConstant(Const('n_IndT', 1));
neuron.AddConstant(Const('n_ActC', 2));
neuron.AddConstant(Const('n_rep1', 5));
neuron.AddConstant(Const('n_ActE', 2));
neuron.AddConstant(Const('n_rep2', 5));
% Production rates
neuron.AddConstant(Const('k_seqrep1', 100));
neuron.AddConstant(Const('k_prodC', 10));
neuron.AddConstant(Const('k_prodD', 10));
neuron.AddConstant(Const('k_prodE', 10));

% Define compositors
drep1dt = neuron.AddCompositor("rep1", in);
dir1dt = neuron.AddCompositor("ir1", 0);
dActCdt = neuron.AddCompositor("ActC", 0);
dActEdt = neuron.AddCompositor("ActE", in * ratio);
drep2dt = neuron.AddCompositor("rep2", 0);
doutdt = neuron.AddCompositor("out", 0);

dIndTdt = neuron.AddCompositor("IndT", threshold);

% Inputs (assumed to be constant and controlled)
neuron.AddPart(Part('rep1', [drep1dt], [Rate('0')]));
neuron.AddPart(Part('ActE', [dActEdt], [Rate('0')]));

% Dynamical equations

neuron.AddPart(Part('out', [doutdt], [Rate('k_prodE * ((K_rep2 ^ n_rep2)/((K_rep2 ^ n_rep2) + rep2 ^ n_rep2)) * (ActE ^ n_ActE / ((K_ActE ^ n_ActE) + ActE ^ n_ActE)) - k_degout * out')]));

neuron.AddPart(Part('rep2', [drep2dt], [Rate('k_prodC * ((K_rep1 ^ n_rep1)/((K_rep1 ^ n_rep1) + rep1 ^ n_rep1)) * (ActC ^ n_ActC / ((K_ActC ^ n_ActC) + ActC ^ n_ActC)) + k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degrep2 * rep2')]));

neuron.AddPart(Part('ir1', [dir1dt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degir1 * ir1')]));

neuron.AddPart(Part('ActC', [dActCdt], [Rate('k_prodD * (IndT ^ n_IndT / ((K_IndT ^ n_IndT) + IndT ^ n_IndT)) - k_degActC * ActC')]));


%% Neuron with the inputs stripped away. Plot the transfer function only

figure()
ax = gca();
hold(ax, 'on')
ax.ColorOrder = jet(20);

for t = 1:20
    s = 60;
    out = zeros(s);
    x_axis = zeros(s);
    for x = 1:s
        ratio = 0.3
        neuron.ChangeInitialValue("IndT", (t-1)/20);
        neuron.ChangeInitialValue("rep1", x);
        neuron.ChangeInitialValue("ActE", x * ratio);
        [T, Y] = neuron.run([0 30]);
        out(x) = Y(15, neuron.CompositorIndex('out'));
        x_axis(x) = x/6;
    end
    ax.ColorOrderIndex = t
    fprintf('index is %d\n',ax.ColorOrderIndex);
    plot(ax, x_axis, out);
end

set(gca, 'FontName', 'Myriad Pro')
title('Transfer function sweep of threshold value, with ratio at 0.3', 'fontsize', 12)


