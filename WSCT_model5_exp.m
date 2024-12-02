%% WISCONSIN SORTING CARD TASK (modifed)

clear all

rng('default')

% Number of time steps per trial
T = 3; 
pObvRule = 3;
pExclRule = 1.2;
% Source deck
deck{1} = {[0, 0, 0, 1]', [0, 0, 0, 1]', [0, 0, 0, 1]'};    % star, yellow, 4
deck{2} = {[0, 0, 1, 0]', [1, 0, 0, 0]', [1, 0, 0, 0]'};    % rectangle, blue, 1
deck{3} = {[1, 0, 0, 0]', [0, 1, 0, 0]', [0, 1, 0, 0]'};    % cicle, red, 2
deck{4} = {[0, 1, 0, 0]', [0, 0, 1, 0]', [0, 0, 1, 0]'};    % triangle, green, 3
deck{5} = {[0, 0, 1, 0]', [1, 0, 0, 0]', [0, 1, 0, 0]'};    % rectangle, blue, 2
deck{6} = {[1, 0, 0, 0]', [0, 1, 0, 0]', [0, 0, 1, 0]'};    % circle, red, 3
deck{7} = {[0, 1, 0, 0]', [0, 0, 1, 0]', [0, 0, 0, 1]'};    % triangle, green, 4
deck{8} = {[0, 0, 0, 1]', [0, 0, 0, 1]', [1, 0, 0, 0]'};    % star, yellow, 1
deck{9} = {[0, 0, 1, 0]', [0, 0, 1, 0]', [0, 0, 0, 1]'};    % rectangle, green, 4
deck{10} = {[1, 0, 0, 0]', [0, 0, 0, 1]', [1, 0, 0, 0]'};   % circle, yellow, 1
deck{11} = {[0, 1, 0, 0]', [1, 0, 0, 0]', [0, 1, 0, 0]'};   % triangle, blue, 2
deck{12} = {[0, 0, 0, 1]', [0, 1, 0, 0]', [0, 0, 1, 0]'};   % star, red, 3
deck{13} = {[1, 0, 0, 0]', [0, 0, 0, 1]', [0, 1, 0, 0]'};   % circle, yellow, 2
deck{14} = {[0, 1, 0, 0]', [1, 0, 0, 0]', [0, 0, 1, 0]'};   % triangle, blue, 3
deck{15} = {[0, 0, 0, 1]', [0, 1, 0, 0]', [0, 0, 0, 1]'};   % star, red, 4
deck{16} = {[0, 0, 1, 0]', [0, 0, 1, 0]', [1, 0, 0, 0]'};   % rectangle, green, 1
deck{17} = {[0, 0, 0, 1]', [1, 0, 0, 0]', [0, 0, 1, 0]'};   % star, blue, 3
deck{18} = {[0, 0, 1, 0]', [0, 1, 0, 0]', [0, 0, 0, 1]'};   % rectangle, red, 4
deck{19} = {[1, 0, 0, 0]', [0, 0, 1, 0]', [1, 0, 0, 0]'};   % circle, green, 1
deck{20} = {[0, 1, 0, 0]', [0, 0, 0, 1]', [0, 1, 0, 0]'};   % triangle, yellow, 2
deck{21} = {[1, 0, 0, 0]', [0, 0, 1, 0]', [0, 0, 1, 0]'};   % circle, green, 3
deck{22} = {[0, 1, 0, 0]', [0, 0, 1, 0]', [0, 0, 0, 1]'};   % triangle, yellow, 4
deck{23} = {[0, 0, 0, 1]', [1, 0, 0, 0]', [1, 0, 0, 0]'};   % star, blue, 1
deck{24} = {[0, 0, 1, 0]', [0, 1, 0, 0]', [0, 1, 0, 0]'};   % rectangle, red, 2

% Target cards
% Card 1
Card{1}{1} = [0, 1, 0, 0]';       % shape = triangle
Card{1}{2} = [0, 1, 0, 0]';       % color = red
Card{1}{3} = [1, 0, 0, 0]';       % number = 1
% Card 2
Card{2}{1} = [0, 0, 0, 1]';       % shape = star
Card{2}{2} = [0, 0, 1, 0]';       % color = green
Card{2}{3} = [0, 1, 0, 0]';       % number = 2
% Card 3
Card{3}{1} = [0, 0, 1, 0]';       % shape = rectangle
Card{3}{2} = [0, 0, 0, 1]';       % color = yellow
Card{3}{3} = [0, 0, 1, 0]';       % number = 3
% Card 4
Card{4}{1} = [1, 0, 0, 0]';       % shape = circle
Card{4}{2} = [1, 0, 0, 0]';       % color = blue
Card{4}{3} = [0, 0, 0, 1]';       % number = 4


% HIDDEN STATE FACTORS
% =========================================================
% Number of hidden state factors
Nf = 6;
% Number of states per factor
Ns(1) = 4; % shape = {circle, triangle, square, star}
Ns(2) = 4; % color = {blue, red, green, yellow}
Ns(3) = 4; % number = {1, 2, 3, 4}
Ns(4) = 4; % rule = {shape, color, number, exclusion}
Ns(5) = 3; % task sequence = {viewing, response, feedback}
Ns(6) = 5; % choice = {card 1, card 2, card 3, card 4, wait}
% Prior for initial states in generative process
D{1} = [0, 0, 0, 1]';         % shape = star
D{2} = [0, 0, 0, 1]';         % color = yellow
D{3} = [0, 0, 0, 1]';         % number = 4
D{4} = [0, 0, 0, 1]';         % rule = exclusion
D{5} = [1, 0, 0]';            % sequence = viewing
D{6} = [0, 0, 0, 0, 1]';      % choice = wait
% Prior for initial states in generative model
d{1} = [0.25, 0.25, 0.25, 0.25]';
d{2} = [0.25, 0.25, 0.25, 0.25]';
d{3} = [0.25, 0.25, 0.25, 0.25]';
d{4} = [pObvRule, pObvRule, pObvRule, pExclRule]';         % lower prior for the exclusion rule
d{5} = [1, 0, 0]';
d{6} = [0, 0, 0, 0, 1]';

% TRANSITION MATRICES
% =========================================================
Nu = 5;                                 % actions = {card 1, card 2, card 3, card 4, wait}
for f=1:4
    B{f}(:, :) = eye(Ns(f));            % features and rule do not change within a trial
end
B{5} = [0 0 1;                          % viewing -> response -> feedback
        1 0 0;
        0 1 0];
B{Nf} = zeros(Ns(Nf), Ns(Nf), Nu);
for u=1:Nu
    B{Nf}(u, :, u) = 1;
end

% Policies
for f=1:Nf-1
    V(:, :, f) = [1 1 1 1;
                  1 1 1 1];                 % features, rule and sequence not controllable
end
V(:, :, Nf) = [5 5 5 5
               1 2 3 4];                    % choice state controllable

% OUTCOME MODALITIES
% =========================================================
% Number of outcome factors
Ng = 5;
% Number of outcomes per factor
No(1) = 4; % shape = {circle, triangle, square, star}
No(2) = 4; % color = {blue, red, green, yellow}
No(3) = 4; % number = {1, 2, 3, 4}
No(4) = 5; % choice = {card 1, card 2, card 3, card 4, wait}
No(5) = 3; % feedback = {incorrect, correct, undecided}

% LIKELIHOOD MAPPING
% ---------------------------------------------------------
for g=1:Ng
    A{g}=zeros([No(g), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5), Ns(6)]);
end
% Shape
A{1}(1, 1, :, :, :, :, :) = 1;     % circle
A{1}(2, 2, :, :, :, :, :) = 1;     % triangle
A{1}(3, 3, :, :, :, :, :) = 1;     % square
A{1}(4, 4, :, :, :, :, :) = 1;     % star
% Color
A{2}(1, :, 1, :, :, :, :) = 1;     % blue
A{2}(2, :, 2, :, :, :, :) = 1;     % red
A{2}(3, :, 3, :, :, :, :) = 1;     % green
A{2}(4, :, 4, :, :, :, :) = 1;     % yellow
% Number
A{3}(1, :, :, 1, :, :, :) = 1;     % 1
A{3}(2, :, :, 2, :, :, :) = 1;     % 2
A{3}(3, :, :, 3, :, :, :) = 1;     % 3
A{3}(4, :, :, 4, :, :, :) = 1;     % 4
% Choice observation
for o=1:No(4)
    A{4}(o, :, :, :, :, :, o) = 1;
end
% Feedback
seq = 3;                           % feedback sequence
rule = 1;                          % shape matching rule
feature = 1;
for shape=1:No(1)
    for choice=1:Ns(6)-1
        if any(shape ~= find(Card{choice}{feature}==1))
            feedback = 1;          % incorrect feedback
            A{5}(feedback, shape, :, :, rule, seq, choice) = 1;

        else
            feedback = 2;          % correct feedback
            A{5}(feedback, shape, :, :, rule, seq, choice) = 1;
        end
    end
end
rule = 2;                       % color matching rule
feature = 2;
for color=1:No(2)
    for choice=1:Ns(6)-1
        if any(color ~= find(Card{choice}{feature}==1))
            feedback = 1;
            A{5}(feedback, :, color, :, rule, seq, choice) = 1;

        else
            feedback = 2;
            A{5}(feedback, :, color, :, rule, seq, choice) = 1;
        end
    end
end
rule = 3;                       % number matching rule
feature = 3;
for num=1:No(3)
    for choice=1:Ns(6)-1
        if any(num ~= find(Card{choice}{feature}==1))
            feedback = 1;
            A{5}(feedback, :, :, num, rule, seq, choice) = 1;

        else
            feedback = 2;
            A{5}(feedback, :, :, num, rule, seq, choice) = 1;
        end
    end
end
rule = 4;                       % no matching feature rule
for choice=1:Ns(6)-1
    for shape=1:No(1)
        for color=1:No(2)
            for num=1:No(3)
                if (shape ~= find(Card{choice}{1})) && (color ~= find(Card{choice}{2})) && (num ~= find(Card{choice}{3}))
                    feedback = 2;
                    A{5}(feedback, shape, color, num, rule, seq, choice) = 1;     
                else
                    feedback = 1;
                    A{5}(feedback, shape, color, num, rule, seq, choice) = 1;
                end
            end
        end
    end
end
% Undecided feedback
A{5}(3, :, :, :, :, :, 5) = 1;              % wait choice
A{5}(3, :, :, :, :, 1, :) = 1;              % viewing sequence
A{5}(3, :, :, :, :, 2, :) = 1;              % response sequence

% PREFERRED OUTCOMES
% ------------------------------------------
loss = 1;
reward = 5;
for i=1:Ng
    C{i} = zeros(No(i), T);
end
C{Ng}(:, T) = [-loss; % Incorrect
                reward;  % Correct
                0]; % Undecided

% Additional parameters
% ---------------------------------------------
% learning rate
eta = 0.5;
% forgetting rate
omega = 0.0;
% expected precision of expected free energy G over policies
beta = 1.0;
% inverse temperature
alpha = 16;

% MDP STRUCTURE
mdp.T = T;                    % Number of time steps
%mdp.U = U;                    % allowable (shallow) policies
mdp.V = V;                    % allowable (deep) policies
mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states

% mdp.a = a; mdp.a_0 = mdp.a;
mdp.d = d; mdp.d0 = mdp.d; mdp.d_0 = mdp.d0;   % enable learning priors over initial states
mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies
mdp.loss = loss;
mdp.reward = reward;
mdp.pObvRule = pObvRule;
mdp.pExclRule = pExclRule;


% We can add labels to states, outcomes, and actions for subsequent plotting:
label.factor{1}   = 'shape';
label.name{1}    = {'circle', 'triangle', 'square', 'star'};
label.factor{2}   = 'color';
label.name{2}    = {'blue', 'red', 'green', 'yellow'};
label.factor{3}   = 'number';
label.name{3}    = {'1', '2', '3', '4'};
label.factor{4}   = 'rule';     
label.name{4}    = {'shape', 'color', 'number', 'exclusion'};
label.factor{5}   = 'sequence';
label.name{5}    = {'viewing', 'response','feedback'};
label.factor{6}   = 'choice';
label.name{6}    = {'card 1','card 2', 'card 3', 'card 4', 'wait'};

label.modality{1}   = 'shape';
label.outcome{1}    = {'circle', 'triangle', 'square', 'star'};
label.modality{2}   = 'color';
label.outcome{2}    = {'blue', 'red', 'green', 'yellow'};
label.modality{3}   = 'number';
label.outcome{3}    = {'1', '2', '3', '4'};
label.modality{4}   = 'choice';
label.outcome{4}    = {'card 1','card 2', 'card 3', 'card 4'};
label.modality{5}   = 'feedback';
label.outcome{5}    = {'incorrect', 'correct', 'undecided'};
for i = 1:Nf
    label.action{i} = {'card1', 'card2', 'card3', 'card4', 'wait'};
end
mdp.label = label;

% MODEL INVERSION
%================================================================
% Participant infos
participant = 381;
N = 99;
isBMR = 'BMR';
% Model reduction
BMR.f = 4;
BMR.pcount = 8;
BMR.thres = 1;
BMR.eps = 0.1;
BMR.rF = {};
BMR.jmin = [];
BMR.rFmin = [];
BMR.rD = [];
BMR.sD = [];
BMR.trials = [];
%BMR1.trial = 30;
OPTIONS.BMR0 = BMR;
% OPTIONS = {};

SIM_TYPE = 'solving';


if strcmp(SIM_TYPE, 'fitting')

    % Get outcomes and actions from experimental subject
    [MDP(1:N)] = deal(mdp);
    Exp_MDP = WSCT_get_data_381(MDP);

    % Get outcomes and actions from simulated subject
    [Sim_Exp(1:N)] = deal(mdp);
    [Sim_Exp, ~] = WSCT_X_tutorial(Sim_Exp, OPTIONS);

    % Free parameters to fit
    params = {'pRuleObv', 'pRuleExcl', 'pcount', 'thres'};

    % Model inversion
    DCM = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, params);

elseif strcmp(SIM_TYPE, 'recovering')

    % Get outcomes and actions from simulated subject
    [Sim_Exp(1:N)] = deal(mdp);
    [Sim_Exp, ~] = WSCT_X_tutorial(Sim_Exp, OPTIONS);

    % Free parameters to fit
    params = {'pRuleObv', 'pRuleExcl', 'pcount', 'thres'};

    % Model inversion
    DCM = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, params);

elseif strcmp(SIM_TYPE, 'solving')
    % Get outcomes and actions from experimental subject
    [MDP(1:N)] = deal(mdp);
    Exp_MDP = WSCT_get_data_381(MDP);

    % Add source cards from experiment to MDP
    MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
    [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);

    % Plotting parameters
    f_act = Nf;
    f_state = 4;
    mod_out = Ng;
    timestep_to_plot = 2;

    % Plot
    WSCT_plot(MDP, f_act, f_state, mod_out, timestep_to_plot);

else
    disp('ERROR: Launch not specified')
end