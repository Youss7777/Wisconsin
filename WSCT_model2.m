%% WISCONSIN SORTING CARD TASK (modifed)

clear

% Number of time steps
T = 2;
% Number of hidden state factors
Nf = 10;
% Number of states per factor
for f=1:Nf-1
    Ns(f) = 2;
end
Ns(Nf) = 3;
% Number of actions
Nu = 2;
% Number of outcome factors
Ng = 1;
% Number of outcomes per factor
No(1) = 3;


% priors for initial states in generative process
% card 1
D{1} = [1, 0]'; % shape of card 1 is a circle
D{2} = [1, 0]'; % color of card 1 is red
D{3} = [1, 0]'; % card 1 has 1 symbol
% card 2
D{4} = [0, 1]'; % shape of card 2 is a triangle
D{5} = [0, 1]'; % color of card 2 is blue
D{6} = [0, 1]'; % card 2 has 2 symbols
% card 3
D{7} = [1, 0]';  % shape of card 3 is a circle
D{8} = [0, 1]'; % color of card 3 is blue
D{9} = [0, 1]'; % card 3 has 2 symbols
% choice state
D{10} = [0, 0, 1]'; % undecided

% priors for initial states in generative model
% card 1
d{1} = [0.25, 0.25]';
d{2} = [0.25, 0.25]';
d{3} = [0.25, 0.25]';
% card 2
d{4} = [0.25, 0.25]';
d{5} = [0.25, 0.25]';
d{6} = [0.25, 0.25]';
% card 3
d{7} = [0.25, 0.25]';
d{8} = [0.25, 0.25]';
d{9} = [0.25, 0.25]';
% choice state
d{10} = [0, 0, 1]'; % agent knows it's undecided


% LIKELIHOOD MAPPING
% --------------------------------------------------------
dims_A = [No(1), Ns(1), Ns(2), Ns(3), Ns(4), Ns(5), Ns(6), Ns(7), Ns(8), Ns(9), Ns(10)];
A{1} = zeros(dims_A);
% outcome = incorrect when shapes don't match
% card 1 is chosen but shapes don't match
A{1}(1, 2, :, :, :, :, :, 1, :, :, 1) = 1;
A{1}(1, 1, :, :, :, :, :, 2, :, :, 1) = 1;
% card 2 is chosen but shapes don't match
A{1}(1, :, :, :, 2, :, :, 1, :, :, 2) = 1;
A{1}(1, :, :, :, 1, :, :, 2, :, :, 2) = 1;

% outcome = correct when shapes match
% choice = card 1 and 1 & 3 shapes match
A{1}(2, 2, :, :, :, :, :, 2, :, :, 1) = 1;
A{1}(2, 1, :, :, :, :, :, 1, :, :, 1) = 1;
% cards 2 is chosen and 2 & 3 shapes match
A{1}(2, :, :, :, 2, :, :, 2, :, :, 2) = 1;
A{1}(2, :, :, :, 1, :, :, 1, :, :, 2) = 1;

% outcome = null when choice_state is undecided
A{1}(3, :, :, :, :, :, :, :, :, :, 3) = 1;

% likelihood mapping in the generative model
a{1} = A{1}*200;
a{1}(1, :, :, :, :, :, :, :, :, :, 1) = 0.25;
a{1}(1, :, :, :, :, :, :, :, :, :, 2) = 0.25;
a{1}(2, :, :, :, :, :, :, :, :, :, 1) = 0.25;
a{1}(2, :, :, :, :, :, :, :, :, :, 2) = 0.25;



% POLICIES (shallow)
% ------------------------------------------------------
% U = zeros(1, 3, 10);
% V = ones(T-1,Nu,Nf);
% V(:, :, 10) = [1 2];
% card features not controllable
for f=1:Nf-1
    U(:, :, f) = [1 1];
end
U(:, :, Nf) = [1 2];
% choice state controllable by actions card1 and card2
% % card 1 features not controllable
% U(:,:,1) = [1 1];
% U(:,:,2) = [1 1];
% U(:,:,3) = [1 1];
% % card 2 features not controllable
% U(:,:,4) = [1 1 1];
% U(:,:,5) = [1 1 1];
% U(:,:,6) = [1 1 1];
% % card 3 features not controllable
% U(:,:,7) = [1 1 1];
% U(:,:,8) = [1 1 1];
% U(:,:,9) = [1 1 1];
% % choice state controllable by actions card1 and card2
% U(:,:,10) = [1 2 2];
% % equivalent way of defining U, per actions and not state factors
% U(:, 1, :) = [10 10 10 10 10 10 10 10 10 10]';
% U(:, 2, :) = [10 10 10 10 10 10 10 10 10 10]';

% TRANSITION MATRICES
% ------------------------------------------------
% card features do not change within a trial
for f=1:Nf-1
    B{f}(:, :) = eye(Ns(f));
end
% only last state factor (choice state) is controllable
B{Nf} = zeros(Ns(Nf), Ns(Nf), Nu);
% choice state
B{10}(1, 1, 1) = 1;
B{10}(1, 2, 1) = 1;
B{10}(1, 3, 1) = 1;
B{10}(2, 1, 2) = 1;
B{10}(2, 2, 2) = 1;
B{10}(2, 3, 2) = 1;



% PREFERRED OUTCOMES
% ------------------------------------------
la = 1;
rs = 4;
C{1}(:,:) =    [0  -la; % Incorrect
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
label.factor{1}   = 'card 1 shape';   label.name{1}    = {'circle','triangle'};
label.factor{2}   = 'card 1 color';     label.name{2}    = {'red', 'blue'};
label.factor{3}   = 'card 1 num';   label.name{3}    = {'1','2'};

label.factor{4}   = 'card 2 shape';   label.name{4}    = {'circle','triangle'};
label.factor{5}   = 'card 2 color';     label.name{5}    = {'red', 'blue'};
label.factor{6}   = 'card 2 num';   label.name{6}    = {'1','2'};

label.factor{7}   = 'card 3 shape';   label.name{7}    = {'circle','triangle'};
label.factor{8}   = 'card 3 color';     label.name{8}    = {'red', 'blue'};
label.factor{9}   = 'card 3 num';   label.name{9}    = {'1','2'};

label.factor{10}   = 'choice state';   label.name{10}    = {'card1', 'card2', 'undecided'};

label.modality{1}   = 'feedback';     label.outcome{1}    = {'incorrect', 'correct', 'null'};

for i = 1:10
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

% trial 2
% MDP(2).D{7} = [0, 1]'; % triangle
% MDP(2).D{8} = [0, 1]'; % blue
% MDP(2).D{9} = [1, 0]'; % 1
% % trial 3
% MDP(3).D{7} = [0, 1]'; % triangle
% MDP(3).D{8} = [1, 0]'; % red
% MDP(3).D{9} = [0, 1]'; % 2
% % trial 4
% MDP(4).D{7} = [1, 0]'; % circle
% MDP(4).D{8} = [1, 0]'; % red
% MDP(4).D{9} = [1, 0]'; % 1


MDP = spm_MDP_VB_X_tutorial(MDP);


WSCT_plot(MDP);




