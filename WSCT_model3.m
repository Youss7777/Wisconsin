%% WISCONSIN SORTING CARD TASK (modifed)

clear all

% Number of time steps per trial
T = 2;
% Card 1
Card{1}{1} = [0, 1]';       % shape = triangle
Card{1}{2} = [0, 1]';       % color = red
Card{1}{3} = [1, 0]';       % number = 1
% Card 2
Card{2}{1} = [1, 0]';       % shape = circle
Card{2}{2} = [1, 0]';       % color = blue
Card{2}{3} = [0, 1]';       % number = 2

% HIDDEN STATE FACTORS
% =========================================================
% Number of hidden state factors
Nf = 5;
% Number of states per factor
Ns(1) = 2; % shape = {circle, triangle}
Ns(2) = 2; % color = {blue, red}
Ns(3) = 2; % number = {1, 2}
Ns(4) = 4; % rule = {shape, color, number, exclusion}
Ns(5) = 3; % choice = {card 1, card 2, undecided}
% Prior for initial states in generative process
D{1} = [0, 1]';         % shape = triangle
D{2} = [1, 0]';         % color = blue
D{3} = [1, 0]';         % number = 1
D{4} = [0, 0, 0, 1]';   % rule = shape
D{5} = [0, 0, 1]';      % choice = undecided
% Prior for initial states in generative model
d{1} = [0.25, 0.25]';
d{2} = [0.25, 0.25]';
d{3} = [0.25, 0.25]';
d{4} = [0.25, 0.25, 0.25, 0.25]';
d{5} = [0, 0, 1]';

% TRANSITION MATRICES
% ----------------------------------------------------------
Nu = 2;                                 % number of actions
for f=1:Nf-1
    B{f}(:, :) = eye(Ns(f));            % features and rule do not change within a trial
end
B{Nf} = zeros(Ns(Nf), Ns(Nf), Nu);
B{Nf}(1, 1, 1) = 1;
B{Nf}(1, 2, 1) = 1;
B{Nf}(1, 3, 1) = 1;
B{Nf}(2, 1, 2) = 1;
B{Nf}(2, 2, 2) = 1;
B{Nf}(2, 3, 2) = 1;
% Actions (shallow policies)
for f=1:Nf-1
    U(:, :, f) = [1 1];                 % features and rule not controllable
end
U(:, :, Nf) = [1 2];                    % choice state controllable

% OUTCOME MODALITIES
% =========================================================
% Number of outcome factors
Ng = 4;
% Number of outcomes per factor
No(1) = 2; % shape = {circle, triangle}
No(2) = 2; % color = {blue, red}
No(3) = 2; % number = {1, 2}
No(4) = 3; % feedback = {incorrect, correct, undecided}

% LIKELIHOOD MAPPING
% ---------------------------------------------------------
A{1}=zeros([No(1), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5)]);
A{2}=zeros([No(2), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5)]);
A{3}=zeros([No(3), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5)]);
A{4}=zeros([No(4), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5)]);
% Shapes
A{1}(1, 1, :, :, :, :) = 1;     % circle
A{1}(2, 2, :, :, :, :) = 1;     % triangle
% Color
A{2}(1, :, 1, :, :, :) = 1;     % blue
A{2}(2, :, 2, :, :, :) = 1;     % red
% Number
A{3}(1, :, :, 1, :, :) = 1;     % 1
A{3}(2, :, :, 2, :, :) = 1;     % 2
% FEEDBACK
rule = 1;               % shape matching rule
feature = 1;
for shape=1:2
    for choice=1:Ns(5)-1
        if any(shape ~= find(Card{choice}{feature}==1))
            feedback = 1;
            A{4}(feedback, shape, :, :, rule, choice) = 1;

        else
            feedback = 2;
            A{4}(feedback, shape, :, :, rule, choice) = 1;
        end
    end
end
rule = 2;               % color matching rule
feature = 2;
for color=1:2
    for choice=1:Ns(5)-1
        if any(color ~= find(Card{choice}{feature}==1))
            feedback = 1;
            A{4}(feedback, :, color, :, rule, choice) = 1;

        else
            feedback = 2;
            A{4}(feedback, :, color, :, rule, choice) = 1;
        end
    end
end
rule = 3;               % number matching rule
feature = 3;
for num=1:2
    for choice=1:Ns(5)-1
        if any(num ~= find(Card{choice}{feature}==1))
            feedback = 1;
            A{4}(feedback, :, :, num, rule, choice) = 1;

        else
            feedback = 2;
            A{4}(feedback, :, :, num, rule, choice) = 1;
        end
    end
end
rule = 4;               % no matching feature rule
for choice=1:Ns(5)-1
    for shape=1:2
        for color=1:2
            for num=1:2
                if any(shape ~= find(Card{choice}{1})) && any(color ~= find(Card{choice}{2})) && any(num ~= find(Card{choice}{3}))
                    feedback = 2;
                    A{4}(feedback, shape, color, num, rule, choice) = 1;     
                else
                    feedback = 1;
                    A{4}(feedback, shape, color, num, rule, choice) = 1;
                end
            end
        end
    end
end
% Null feedback when undecided choice
A{4}(3, :, :, :, :, 3) = 1;


% Likelihood mapping in the generative model
for i=1:Ng
    a{i}=A{i}*200;
end
a{4}(1, :, :, :, 1) = 0.25;
a{4}(2, :, :, :, 1) = 0.25;
a{4}(1, :, :, :, 2) = 0.25;
a{4}(2, :, :, :, 2) = 0.25;

% PREFERRED OUTCOMES
% ------------------------------------------
la = 1;
rs = 4;
for i=1:Ng-1
    C{i}(:,:) = [0 0]';
end
C{4}(:,:) =    [0  -la; % Incorrect
                0 rs;  % Correct
                0  0]; % Null

% Additional parameters
% ---------------------------------------------
% learning rate
eta = 0.5;
% forgetting rate
omega = 0.0;
% expected precision of expected free energy G over policies
beta = 1.0;
% inverse temperature
alpha = 512; % deterministic action (always the most probable)

% MDP STRUCTURE
mdp.T = T;                    % Number of time steps
mdp.U = U;                    % allowable (shallow) policies
% mdp.V = V;                    % allowable (deep) policies
mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states

mdp.d = d; mdp.d_0 = mdp.d;   % enable learning priors over initial states
mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies

mdp.a = a; mdp.a_0 = mdp.a;
% mdp.s = s;
% mdp.o = o;

% We can add labels to states, outcomes, and actions for subsequent plotting:
label.factor{1}   = 'shape';
label.name{1}    = {'circle', 'triangle'};
label.factor{2}   = 'color';
label.name{2}    = {'blue', 'red'};
label.factor{3}   = 'number';
label.name{3}    = {'1', '2'};
label.factor{4}   = 'rule';     
label.name{4}    = {'shape', 'color', 'number', 'exclusion'};
label.factor{5}   = 'choice';
label.name{5}    = {'card 1','card 2', 'undecided'};

label.modality{1}   = 'shape';
label.outcome{1}    = {'circle', 'triangle'};
label.modality{2}   = 'color';
label.outcome{2}    = {'blue', 'red'};
label.modality{3}   = 'number';
label.outcome{3}    = {'1', '2'};
label.modality{4}   = 'feedback';
label.outcome{4}    = {'incorrect', 'correct', 'null'};

for i = 1:Nf
    label.action{i} = {'card1', 'card2'};
end
mdp.label = label;

clear beta
clear alpha
clear eta
clear omega
clear la
clear rs % We clear these so we can re-specify them in later simulations



% Multiple trials simulation
N = 40; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

% Changing features
for i=2:N
    for feature=1:3
        MDP(i).D{feature} = zeros(Ns(2), 1);
        rand_idx = randi([1, 2]);
        MDP(i).D{feature}(rand_idx) = 1;
    end
end


MDP = spm_MDP_VB_X_tutorial(MDP);

f_act = Nf;
f_state = 4;
mod_out = Ng;
timestep_of_trial_to_plot = 1;

WSCT_plot(MDP, f_act, f_state, mod_out, timestep_of_trial_to_plot);




