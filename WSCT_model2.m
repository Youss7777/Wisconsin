%% WISCONSIN SORTING CARD TASK (modifed)

clear all

% Number of time steps per trial
T = 2;

% HIDDEN STATE FACTORS
% =========================================================
% Number of hidden state factors
Nf = 3;
% We assume: Card 1 = {triangle x red x 1}, Card 2 = {circle x blue x 2}
% Number of states per factor
Ns(1) = 4; % rule = {shape, color, number, exclusion}
Ns(2) = 8; % conjunction = {circle x blue x 1, circle x blue x 2,
           %   (card 3)    circle x red x 1, circle x red x 2,
           %               triangle x blue x 1, triangle x blue x 2,
           %               triangle x red x 1, triangle x red x 2}
Ns(3) = 3; % choice = {card 1, card 2, undecided}
% Prior for initial states in generative process
D{1} = [0, 0, 0, 1]';                     % rule = shape
D{2} = [0, 0, 0, 0, 1, 0, 0, 0]';         % card 1 = triangle x blue x 1
D{3} = [0, 0, 1]';                        % choice = undecided
% Prior for initial states in generative model
d{1} = [0.25, 0.25, 0.25, 0.25]';
d{2} = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]';
d{3} = [0, 0, 1]';
% TRANSITION MATRICES
% --------------------------------------------------------
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
% --------------------------------------------------------
A{1}=zeros([No(1), Ns(1), Ns(2), Ns(3)]);
A{2}=zeros([No(2), Ns(1), Ns(2), Ns(3)]);
A{3}=zeros([No(3), Ns(1), Ns(2), Ns(3)]);
A{4}=zeros([No(4), Ns(1), Ns(2), Ns(3)]);

% Shape = circle
A{1}(1, :, 1, :) = 1;
A{1}(1, :, 2, :) = 1;
A{1}(1, :, 3, :) = 1;
A{1}(1, :, 4, :) = 1;
% Shape = triangle
A{1}(2, :, 5, :) = 1;
A{1}(2, :, 6, :) = 1;
A{1}(2, :, 7, :) = 1;
A{1}(2, :, 8, :) = 1;
% Color = blue
A{2}(1, :, 1, :) = 1;
A{2}(1, :, 2, :) = 1;
A{2}(1, :, 5, :) = 1;
A{2}(1, :, 6, :) = 1;
% Color = red
A{2}(2, :, 3, :) = 1;
A{2}(2, :, 4, :) = 1;
A{2}(2, :, 7, :) = 1;
A{2}(2, :, 8, :) = 1;
% Num = 1
A{3}(1, :, 1, :) = 1;
A{3}(1, :, 3, :) = 1;
A{3}(1, :, 5, :) = 1;
A{3}(1, :, 7, :) = 1;
% Num 2
A{3}(2, :, 2, :) = 1;
A{3}(2, :, 4, :) = 1;
A{3}(2, :, 6, :) = 1;
A{3}(2, :, 8, :) = 1;
% Feedback = null
A{4}(3, :, :, 3) = 1;        % choice = undecided
% RULE = MATCHING BY SHAPE
% feedback = incorrect
A{4}(1, 1, 1, 1) = 1;        % card 1 chosen but shape = circle
A{4}(1, 1, 2, 1) = 1;
A{4}(1, 1, 3, 1) = 1;
A{4}(1, 1, 4, 1) = 1;
A{4}(1, 1, 5, 2) = 1;        % card 2 chosen but shape = triangle
A{4}(1, 1, 6, 2) = 1;
A{4}(1, 1, 7, 2) = 1;
A{4}(1, 1, 8, 2) = 1;
% feedback = correct
A{4}(2, 1, 5, 1) = 1;        % card 1 chosen and shape = triangle
A{4}(2, 1, 6, 1) = 1;
A{4}(2, 1, 7, 1) = 1;
A{4}(2, 1, 8, 1) = 1;
A{4}(2, 1, 1, 2) = 1;        % card 2 chosen and shape = cicle
A{4}(2, 1, 2, 2) = 1;
A{4}(2, 1, 3, 2) = 1;
A{4}(2, 1, 4, 2) = 1;
% RULE = MATCHING BY COLOR
% feedback = incorrect
A{4}(1, 2, 1, 1) = 1;        % card 1 chosen but color = blue
A{4}(1, 2, 2, 1) = 1;
A{4}(1, 2, 5, 1) = 1;
A{4}(1, 2, 6, 1) = 1;
A{4}(1, 2, 3, 2) = 1;        % card 2 chosen but color = red
A{4}(1, 2, 4, 2) = 1;
A{4}(1, 2, 7, 2) = 1;
A{4}(1, 2, 8, 2) = 1;
% feedback = correct
A{4}(2, 2, 3, 1) = 1;        % card 1 chosen and color = red
A{4}(2, 2, 4, 1) = 1;
A{4}(2, 2, 7, 1) = 1;
A{4}(2, 2, 8, 1) = 1;
A{4}(2, 2, 1, 2) = 1;        % card 2 chosen and color = blue
A{4}(2, 2, 2, 2) = 1;
A{4}(2, 2, 5, 2) = 1;
A{4}(2, 2, 6, 2) = 1;
% RULE = MATCHING BY NUMBER
% feedback = incorrect
A{4}(1, 3, 2, 1) = 1;        % card 1 chosen but num = 2
A{4}(1, 3, 4, 1) = 1;
A{4}(1, 3, 6, 1) = 1;
A{4}(1, 3, 8, 1) = 1;
A{4}(1, 3, 1, 2) = 1;        % card 2 chosen but num = 1
A{4}(1, 3, 3, 2) = 1;
A{4}(1, 3, 5, 2) = 1;
A{4}(1, 3, 7, 2) = 1;
% feedback = correct
A{4}(2, 3, 1, 1) = 1;        % card 1 chosen and num = 1
A{4}(2, 3, 3, 1) = 1;
A{4}(2, 3, 5, 1) = 1;
A{4}(2, 3, 7, 1) = 1;
A{4}(2, 3, 2, 2) = 1;        % card 2 chosen and num = 2
A{4}(2, 3, 4, 2) = 1;
A{4}(2, 3, 6, 2) = 1;
A{4}(2, 3, 8, 2) = 1;
% RULE = NO MATCHING FEATURES (exclusion)
% feedback = incorrect
A{4}(1, 4, 1, 1) = 1;        % card 1 chosen but at least one same feature
A{4}(1, 4, 3, 1) = 1;
A{4}(1, 4, 4, 1) = 1;
A{4}(1, 4, 5, 1) = 1;
A{4}(1, 4, 6, 1) = 1;
A{4}(1, 4, 7, 1) = 1;
A{4}(1, 4, 8, 1) = 1;
A{4}(1, 4, 1, 2) = 1;        % card 2 chosen but at least one same  feature
A{4}(1, 4, 2, 2) = 1;
A{4}(1, 4, 3, 2) = 1;
A{4}(1, 4, 4, 2) = 1;
A{4}(1, 4, 5, 2) = 1;
A{4}(1, 4, 6, 2) = 1;
A{4}(1, 4, 8, 2) = 1;
% feedback = correct
A{4}(2, 4, 2, 1) = 1;        % card 1 chosen and shape=circle, color=blue, num=2
A{4}(2, 4, 7, 2) = 1;        % card 2 chosen and shape=triangle, color=red, num=1
% Likelihood mapping in the generative model
for i=1:Ng
    a{i}=A{i}*200;
end
a{4}(1, :, :, 1) = 0.25;
a{4}(2, :, :, 1) = 0.25;
a{4}(1, :, :, 2) = 0.25;
a{4}(2, :, :, 2) = 0.25;


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

% True states and observations at time step 1
% s = [1; 1; 1; 2; 2; 2; 1; 2; 2; 3];
% o = 3;


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
label.factor{1}   = 'rule';     
label.name{1}    = {'shape', 'color', 'number', 'exclusion'};
label.factor{2}   = 'conjunctions of target card';   
label.name{2}    = {'circle x blue x 1, circle x blue x 2';
                    'circle x red x 1, circle x red x 2';
                    'triangle x blue x 1, triangle x blue x 2';
                    'triangle x red x 1, triangle x red x 2'};
label.factor{3}   = 'choice';
label.name{3}    = {'card 1','card 2', 'undecided'};

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
N = 20; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

% Changing features
for i=2:N
    MDP(i).D{2} = zeros(Ns(2), 1);
    rand_idx = randi([1, 8]);
    MDP(i).D{2}(rand_idx) = 1;
end


MDP = spm_MDP_VB_X_tutorial(MDP);

f_act = 3;
f_state = 1;
mod_out = 4;
timestep_of_trial_to_plot = 1;

WSCT_plot(MDP, f_act, f_state, mod_out, timestep_of_trial_to_plot);




