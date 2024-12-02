function [DCM] = WSCT_Estimate_parameters(DCM)

% MDP inversion using Variational Bayes
% FORMAT [DCM] = spm_dcm_mdp(DCM)
%
% Expects:
%--------------------------------------------------------------------------
% DCM.MDP       % MDP structure specifying a generative model
% DCM.OPTIONS   % OPTIONS structure (e.g. for BMR)
% DCM.field     % parameter (field) names to optimise
% DCM.U         % cell array of outcomes (stimuli)
% DCM.Y         % cell array of responses (action)
%
% Returns:
%--------------------------------------------------------------------------
% DCM.M     % generative model (DCM)
% DCM.Ep    % Conditional means (structure)
% DCM.Cp    % Conditional covariances
% DCM.F     % (negative) Free-energy bound on log evidence
% 
% This routine inverts (cell arrays of) trials specified in terms of the
% stimuli or outcomes and subsequent choices or responses. It first
% computes the prior expectations (and covariances) of the free parameters
% specified by DCM.field. These parameters are log scaling parameters that
% are applied to the fields of DCM.MDP. 
%
% If there is no learning implicit in multi-trial games, only unique trials
% (as specified by the stimuli), are used to generate (subjective)
% posteriors over choice or action. Otherwise, all trials are used in the
% order specified. The ensuing posterior probabilities over choices are
% used with the specified choices or actions to evaluate their log
% probability. This is used to optimise the MDP (hyper) parameters in
% DCM.field using variational Laplace (with numerical evaluation of the
% curvature).
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_mdp.m 7120 2017-06-20 11:30:30Z spm $

% OPTIONS
%--------------------------------------------------------------------------
ALL = false;

% Here we specify prior expectations (for parameter means and variances)
%--------------------------------------------------------------------------
prior_variance = 1/4; % smaller values will lead to a greater complexity 
                      % penalty (posteriors will remain closer to priors)
prior_alpha = 16;
prior_beta = 1;
prior_loss = 1;
prior_reward = 5;
prior_eta = 0.5;
prior_omega = 0.5;
prior_pRuleObv = 3;
prior_pRuleExcl = 1.2;
prior_pcount = 8;
prior_thres = 1;

for i = 1:length(DCM.field)
    field = DCM.field{i};
    try
        param = DCM.MDP.(field);
        param = double(~~param);
    catch
        param = 1;
    end
    if ALL
        pE.(field) = zeros(size(param));
        pC{i,i}    = diag(param);
    else
        if strcmp(field,'alpha')
            pE.(field) = log(prior_alpha);          % in log-space (to keep positive)
            pC{i,i}    = prior_variance;
        elseif strcmp(field,'beta')
            pE.(field) = log(prior_beta);           % in log-space (to keep positive)
            pC{i,i}    = prior_variance;
        elseif strcmp(field,'loss')
            pE.(field) = log(prior_loss);           % in log-space (to keep positive)
            pC{i,i}    = prior_variance;
        elseif strcmp(field,'reward')
            pE.(field) = log(prior_reward);           % in log-space (to keep positive)
            pC{i,i}    = prior_variance;
        elseif strcmp(field,'eta')
            pE.(field) = log(prior_eta/(1-prior_eta)); % in logit-space - bounded between 0 and 1
            pC{i,i}    = prior_variance;
        elseif strcmp(field,'omega')
            pE.(field) = log(prior_omega/(1-prior_omega)); % in logit-space - bounded between 0 and 1
            pC{i,i}    = prior_variance;
        elseif strcmp(field, 'pRuleObv')
            pE.(field) = log(prior_pRuleObv);
            pC{i,i}    = prior_variance;
        elseif strcmp(field, 'pRuleExcl')
            pE.(field) = log(prior_pRuleExcl);
            pC{i,i}    = prior_variance;
        elseif strcmp(field, 'pcount')
            pE.(field) = log(prior_pcount);
            pC{i,i}    = prior_variance;
        elseif strcmp(field, 'thres')
            pE.(field) = log(prior_thres);
            pC{i,i}    = prior_variance;
        else
            pE.(field) = 0;                % if it can take any negative or positive value
            pC{i,i}    = prior_variance;
        end
    end
end

pC      = spm_cat(pC);

% model specification
%--------------------------------------------------------------------------
M.L     = @(P,M,U,Y)spm_mdp_L(P,M,U,Y);  % log-likelihood function
M.pE    = pE;                            % prior means (parameters)
M.pC    = pC;                             % prior variance (parameters)
M.mdp   = DCM.mdp;                       % MDP structure
M.OPTIONS = DCM.OPTIONS;                % OPTIONS structure

% Variational Laplace
%--------------------------------------------------------------------------
[Ep,Cp,F] = spm_nlsi_Newton(M,DCM.U,DCM.Y); % This is the actual fitting routine

% Store posterior distributions and log evidence (free energy)
%--------------------------------------------------------------------------
DCM.M   = M;  % Generative model
DCM.Ep  = Ep; % Posterior parameter estimates
DCM.Cp  = Cp; % Posterior variances and covariances
DCM.F   = F;  % Free energy of model fit


return

function L = spm_mdp_L(P,M,U,Y)
% log-likelihood function
% FORMAT L = spm_mdp_L(P,M,U,Y)
% P    - parameter structure
% M    - generative model
% U    - inputs
% Y    - observed repsonses
%
% This function runs the generative model with a given set of parameter
% values, after adding in the observations and actions on each trial
% from (real or simulated) participant data. It then sums the
% (log-)probabilities (log-likelihood) of the participant's actions under the model when it
% includes that set of parameter values. The variational Bayes fitting
% routine above uses this function to find the set of parameter values that maximize
% the probability of the participant's actions under the model (while also
% penalizing models with parameter values that move farther away from prior
% values).
%__________________________________________________________________________

if ~isstruct(P); P = spm_unvec(P,M.pE); end


% Here we re-transform parameter values out of log- or logit-space when 
% inserting them into the model to compute the log-likelihood
%--------------------------------------------------------------------------
mdp   = M.mdp;
OPTIONS = M.OPTIONS;
field = fieldnames(M.pE);
for i = 1:length(field)
    if strcmp(field{i},'alpha')
        mdp.(field{i}) = exp(P.(field{i}));
        disp('alpha: '), disp(mdp.alpha)
    elseif strcmp(field{i},'beta')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'loss')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'reward')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'eta')
        mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        disp('eta: '), disp(mdp.eta)
    elseif strcmp(field{i},'omega')
        mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
    elseif strcmp(field{i}, 'pRuleObv')
        mdp.(field{i}) = exp(P.(field{i}));
        disp('pRuleObv: '), disp(mdp.pRuleObv)
    elseif strcmp(field{i}, 'pRuleExcl')
        mdp.(field{i}) = exp(P.(field{i}));
        disp('pRuleExcl: '), disp(mdp.pRuleExcl)
    elseif strcmp(field{i}, 'pcount')
        OPTIONS.BMR0.(field{i}) = exp(P.(field{i}));
        disp('pcount: '), disp(OPTIONS.BMR0.pcount)
    elseif strcmp(field{i}, 'thres')
        OPTIONS.BMR0.(field{i}) = exp(P.(field{i}));
        disp('thres: '), disp(OPTIONS.BMR0.thres)
    else
        mdp.(field{i}) = exp(P.(field{i}));
    end
end

% place MDP in trial structure
%--------------------------------------------------------------------------

if isfield(M.pE,'loss')&&isfield(M.pE,'reward')
    mdp.C{5}(1, 3) = mdp.loss;
    mdp.C{5}(2, 3) = mdp.reward;
elseif isfield(M.pE,'loss')
    mdp.C{5}(1, 3) = -mdp.loss;
elseif isfield(M.pE,'reward')
    mdp.C{5}(2, 3) = mdp.reward;
end

if isfield(M.pE, 'pRuleObv')
    mdp.d{4}(1:3) = mdp.pRuleObv;
    mdp.d0 = mdp.d; mdp.d_0 = mdp.d0;
end

if isfield(M.pE, 'pRuleExcl')
    mdp.d{4}(4)= mdp.pRuleExcl;
    mdp.d0 = mdp.d; mdp.d_0 = mdp.d0;
end

if isfield(M.pE, 'pcount')
   OPTIONS.BMR0.pcount = OPTIONS.BMR0.pcount;
end
if isfield(M.pE, 'thres')
   OPTIONS.BMR0.thres = OPTIONS.BMR0.thres;
end

j = 1:numel(U); % observations for each trial (U = DCM.U = outcomes)
n = numel(j);   % number of trials

[MDP(1:n)] = deal(mdp);  % Create MDP with number of specified trials
[MDP.o]    = deal(U{j}); % Add observations in each trial
% Add correct source cards from experiment to each trial
MDP = WSCT_draw_from_exp_obs(MDP, U);
controllable_factor = 5;

% solve MDP and accumulate log-likelihood
%--------------------------------------------------------------------------
[MDP, OPTIONS]   = WSCT_X_tutorial(MDP, OPTIONS); % run model with possible parameter values


L     = 0; % start (log) probability of actions given the model at 0

for i = 1:numel(Y) % Get probability of true actions for each trial
    for j = 1:numel(Y{1}(controllable_factor,:)) % Only get probability of the second (controllable) state factor

        L = L + log(MDP(i).P(:,:,:,:,:,Y{i}(controllable_factor,j),j)+ eps); % sum the (log) probabilities of each action
                                                   % given a set of possible parameter values
    end
end 

clear('MDP')

fprintf('LL: %f \n',L)

