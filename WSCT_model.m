%% WISCONSIN SORTING CARD TASK (modifed)
clear all
rng('default')
% Number of time steps per trial
T = 3; 
% PARAMETERS
% ---------------------------------------------
% inverse temperature
alpha = 12;
% learning rate
eta = 0.5;
% forgetting rate
omega = 0.0;
% expected precision of expected free energy G over policies
beta = 1.0;
% loss and reward
loss = 1;
reward = 5;
% prior concentration parameters for rules
pRuleObv = 2;
pRuleExcl = 1;
% BMR
BMR.pcount = 2.5;
BMR.thres = 0.7;
isBMR = 'BMR';
% Choose type of simulation
SIM_TYPE = 'solving';
SOURCE_CARDS = 'random';    % random (from unambiguous source deck), 381, ambiguous (features change randomly)
% Participant infos
participant = 381;
N_exp = 99;
N = 50;
N_agents = 60;
% Free parameters to fit
params = {'alpha', 'eta', 'pRuleExcl'};
FIT_PARAMS.alpha = alpha;
FIT_PARAMS.eta = eta;
FIT_PARAMS.omega = omega;
FIT_PARAMS.beta = beta;
FIT_PARAMS.loss = loss;
FIT_PARAMS.reward = reward;
FIT_PARAMS.pRuleObv = pRuleObv;
FIT_PARAMS.pRuleExcl = pRuleExcl;
FIT_PARAMS.pcount = BMR.pcount;
FIT_PARAMS.thres = BMR.thres;

% SOURCE DECK
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

% TARGET CARDS
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
d{4} = [pRuleObv, pRuleObv, pRuleObv, pRuleExcl]';         % lower prior for the exclusion rule
d{5} = [1, 0, 0]';
d{6} = [0, 0, 0, 0, 1]';

% TRANSITION MATRICES
% =========================================================
Nu = 5;                                 % actions = {card 1, card 2, card 3, card 4, wait}
for f=1:4
    B{f}(:, :) = eye(Ns(f));            % features do not change within a trial
end
% B{4}(:,:) = [0.97 0.01 0.01 0.01;       % adding some stochasticity to model ambiguity of humans' models
%              0.01 0.97 0.01 0.01;
%              0.01 0.01 0.97 0.01;
%              0.01 0.01 0.01 0.97];
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
for i=1:Ng
    C{i} = zeros(No(i), T);
end
C{Ng}(:, T) = [-loss; % Incorrect
                reward;  % Correct
                0]; % Undecided


% MDP STRUCTURE
%==================================================================
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
mdp.pRuleObv = pRuleObv;
mdp.pRuleExcl = pRuleExcl;

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

% Plotting parameters
f_act = Nf;
f_state = 4;
mod_out = Ng;
timestep_to_plot = 2;

% MODEL INVERSION
%==================================================================
% Model reduction
BMR.f = 4;
BMR.eps = 0.1;
for i=1:size(d{4}, 1)
    BMR.rF{i} = [];
end
BMR.jmin = [];
BMR.rFmin = [];
BMR.rD = [];
BMR.sD = [];
BMR.trials = [];
%BMR.trial = 82;
if strcmp(isBMR, 'BMR')
    OPTIONS.BMR0 = BMR;
else
    OPTIONS = {};
end


if strcmp(SIM_TYPE, 'fitting')
    N = N_exp;
    % Get outcomes and actions from experimental subject
    [MDP(1:N)] = deal(mdp);
    if participant == 768
        Exp_MDP = WSCT_get_data_768(MDP);
    elseif participant == 381
        Exp_MDP = WSCT_get_data_381(MDP);        
    elseif participant == 416
        Exp_MDP = WSCT_get_data_416(MDP);        
    elseif participant == 828
        Exp_MDP = WSCT_get_data_828(MDP);
    elseif participant == 3
        Exp_MDP = WSCT_get_data_3(MDP);    
    end
    disp(['participant: ' num2str(participant)])
    % Model inversion
    DCM = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, params, FIT_PARAMS);
    disp('Saving DCM...')
    filename = ['WSCT_DCM_deep_' num2str(participant) '_' isBMR '_ ' SIM_TYPE '.mat'];
    save(filename, 'DCM');
    disp('DCM saved!')

elseif strcmp(SIM_TYPE, 'recovering')
    N = N_exp;
    % Get outcomes and actions from simulated subject
    [Sim_Exp(1:N)] = deal(mdp);
    [Sim_Exp, ~] = WSCT_X_tutorial(Sim_Exp, OPTIONS);
    disp(['participant: ' num2str(participant)])
    % Model inversion
    DCM = WSCT_model_inversion(mdp, OPTIONS, Sim_Exp, params, FIT_PARAMS);
    % Saving DCM
    disp('Saving DCM...')
    filename = ['WSCT_DCM_deep_' num2str(participant) '_' isBMR '_' SIM_TYPE '.mat'];
    save(filename, 'DCM');
    disp('DCM saved!')

elseif strcmp(SIM_TYPE, 'comparison')
    N = N_exp;
    [MDP(1:N)] = deal(mdp);
    if participant == 768
        Exp_MDP = WSCT_get_data_768(MDP);
    elseif participant == 381
        Exp_MDP = WSCT_get_data_381(MDP);        
    elseif participant == 416
        Exp_MDP = WSCT_get_data_416(MDP);        
    elseif participant == 828
        Exp_MDP = WSCT_get_data_828(MDP);
    elseif participant == 3
        Exp_MDP = WSCT_get_data_3(MDP);
    end
    disp(['participant: ' num2str(participant)])
    % BMR
    clear isBMR
    isBMR = 'BMR';
    DCM_BMR = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, params, FIT_PARAMS);
    disp('Saving DCM...')
    filename_BMR = ['WSCT_DCM_deep_' num2str(participant) '_' isBMR '_' SIM_TYPE '.mat'];
    save(filename_BMR, 'DCM_BMR');
    disp('DCM saved!')

    % No BMR
    clear OPTIONS isBMR
    OPTIONS = {};
    isBMR = 'noBMR';
    DCM_noBMR = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, params, FIT_PARAMS);
    disp('Saving DCM...')
    filename_noBMR = ['WSCT_DCM_deep_' num2str(participant) '_' isBMR '_' SIM_TYPE '.mat'];
    save(filename_noBMR, 'DCM_noBMR');
    disp('DCM saved!')

    % Model comparison
    DCM_381_BMR = load(filename_BMR);
    F_381_BMR = DCM_381_BMR.DCM_BMR.F;
    DCM_381_noBMR = load(filename_noBMR);
    F_381_noBMR = DCM_381_noBMR.DCM_noBMR.F;
    [alpha,exp_r,xp,pxp,bor] = spm_BMS([F_381_BMR F_381_noBMR]);
    disp(['pxp= ' mat2str(pxp)])

elseif strcmp(SIM_TYPE, 'performance')
    [MDP(1:N)] = deal(mdp);
    if strcmp(SOURCE_CARDS, 'random')
        Rand_deck_MDP = WSCT_draw_from_deck_random(MDP, deck);
        WSCT_performance(Rand_deck_MDP, OPTIONS, N, N_agents, f_state)
    end


elseif strcmp(SIM_TYPE, 'solving')
    if strcmp(SOURCE_CARDS, '381')
        N = N_exp;
        [MDP(1:N)] = deal(mdp);
        % Add source cards from experiment to MDP
        Exp_MDP = WSCT_get_data_381(MDP);
        disp(['participant: ' num2str(participant)])
        % Get source cards from subject observations
        MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
        [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);
        WSCT_nb_correct_pred_act(MDP, Exp_MDP)        
    elseif strcmp(SOURCE_CARDS, '416')
        N = N_exp;
        [MDP(1:N)] = deal(mdp);
        % Add source cards from experiment to MDP
        Exp_MDP = WSCT_get_data_416(MDP);
        disp(['participant: ' num2str(participant)])
        % Get source cards from subject observations
        MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
        [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);
        WSCT_nb_correct_pred_act(MDP, Exp_MDP)        
    elseif strcmp(SOURCE_CARDS, '768')
        N = N_exp;
        [MDP(1:N)] = deal(mdp);
        % Add source cards from experiment to MDP
        Exp_MDP = WSCT_get_data_768(MDP);
        disp(['participant: ' num2str(participant)])
        % Get source cards from subject observations
        MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
        [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);
        WSCT_nb_correct_pred_act(MDP, Exp_MDP)        
    elseif strcmp(SOURCE_CARDS, '3')
        N = N_exp;
        [MDP(1:N)] = deal(mdp);
        % Add source cards from experiment to MDP
        Exp_MDP = WSCT_get_data_3(MDP);
        disp(['participant: ' num2str(participant)])
        % Get source cards from subject observations
        MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
        [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);
        WSCT_nb_correct_pred_act(MDP, Exp_MDP)        
    elseif strcmp(SOURCE_CARDS, '828')
        N = N_exp;
        [MDP(1:N)] = deal(mdp);
        % Add source cards from experiment to MDP
        Exp_MDP = WSCT_get_data_828(MDP);
        disp(['participant: ' num2str(participant)])
        % Get source cards from subject observations
        MDP = WSCT_draw_from_exp_obs(MDP, {Exp_MDP.o});
        [MDP, OPTIONS] = WSCT_X_tutorial(MDP, OPTIONS);
        WSCT_nb_correct_pred_act(MDP, Exp_MDP)
    elseif strcmp(SOURCE_CARDS, 'random')
        [MDP(1:N)] = deal(mdp);        
        % Add source cards from deck randomly        
        Rand_deck_MDP = WSCT_draw_from_deck_random(MDP, deck);
        [MDP, OPTIONS] = WSCT_X_tutorial(Rand_deck_MDP, OPTIONS);
    elseif strcmp(SOURCE_CARDS, 'ambiguous')
        [MDP(1:N)] = deal(mdp);        
        % Add source cards with random features       
        Rand_MDP = WSCT_draw_random(MDP, Ns);
        [MDP, OPTIONS] = WSCT_X_tutorial(Rand_MDP, OPTIONS);
    else
        disp('Error. Source cards not specified.')
    end


    % Plot
    surprise = WSCT_plot(MDP, OPTIONS, f_act, f_state, mod_out, timestep_to_plot);

else
    disp('ERROR: Launch not specified')
end


function [perc] = WSCT_nb_correct_pred_act(MDP, Exp_MDP)
s = 0;
for i=1:length(MDP)
    ss = sum(Exp_MDP(i).u(6,2)==MDP(i).u(6,2));
    s = s + ss;
end
perc = s/length(MDP);
disp(['Percentage of correct predicted actions: ' num2str(perc)])
end