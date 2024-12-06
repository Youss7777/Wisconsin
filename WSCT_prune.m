function [MDP, BMR] = WSCT_prune(MDP, BMR)
% FORMAT [MDP, BMR] = WSCT_prune(MDP, BMR)
% 
% MDP - Markov Decision Process structure containing priors and posteriors
% BMR - Structure containing parameters for Bayesian model reduction
%
% This function performs Bayesian model reduction on the concentration
% parameters in MDP.d{f}, ensuring that concentration parameters remain
% non-negative and valid.
%__________________________________________________________________________

% Extract parameters from BMR structure
%--------------------------------------------------------------------------
f = BMR.f;    % Hidden factor to integrate over
pcount = BMR.pcount;    % Prior counts added during model reduction
thres = BMR.thres;    % Threshold for Bayesian model reduction
eps = BMR.eps;

% Extract prior and posterior concentration parameters
%--------------------------------------------------------------------------
qD = MDP.d{f};      % Posterior concentration parameters
pD = MDP.d_0{f};    % Prior concentration parameters

% Initialize reduced posterior and prior
%--------------------------------------------------------------------------
sD = qD;
rD = pD;

% Generate model space: additional concentration parameters (i.e., precision)
%--------------------------------------------------------------------------
% n = numel(pD);
% motifs = ff2n(n);           % All combinations of binary vectors of dimension n
% motifs(1, : ) = [];         % Remove all-zero vector
% motifs(end, : ) = [];       % Remove all-ones vector

% motifs(1, :) = [1 0 0 0]';
% motifs(2, :) = [0 1 0 0]';
% motifs(3, :) = [0 0 1 0]';
% motifs(4, :) = [0 0 0 1]';

% Identify indices of non-zero prior parameters
%--------------------------------------------------------------------------
j = find(pD > eps);       % Indices where prior is non-zero
p = pD; %p = pD(j);
q = qD; %q = qD(j);
% Only proceed if number of parameters is greater than 1
%--------------------------------------------------------------------------
if length(j) > 1
    % disp(['pD = ' mat2str(pD)])
    % disp(['qD = ' mat2str(qD)])
    % disp(['j = ' mat2str(j)])
    % disp(['r (p) = ' mat2str(p)])
    F = zeros(length(j), 1);
    for i = 1:length(j)
        r = p;
        r(i) = r(i) + pcount;       % r(i) + pcount normally
        F(i) = spm_MDP_log_evidence(q,p,r);
        BMR.rF{i} = [BMR.rF{i} F(i)];             % store expected free energy reduction
        disp(['rF{' num2str(i) '} = ' mat2str(BMR.rF{i})])

        disp(['p = ' mat2str(p)])
        disp(['q = ' mat2str(q)])
        disp(['r = ' mat2str(r)])
        disp(['F(i) = ' num2str(F(i))])
    end
    [Fmin, imin] = min(F);
    jmin = j(imin);   % Index in pD
    
    if Fmin < -thres
        sD(:) = eps;
        rD(:) = eps;
        sD(jmin) = sum(q);
        rD(jmin) = sum(p);
    
        BMR.applied     = true;
        BMR.jmin    = [BMR.jmin jmin];
        BMR.rFmin   = [BMR.rFmin Fmin];
        BMR.rD      = [BMR.rD rD];
        BMR.sD      = [BMR.sD sD];
    
        disp(['r = ' mat2str(r)])
        disp(['F(i) = ' num2str(F(i))])
        disp('Selected reduced prior (rD):'), disp(rD)
        disp('Updated posterior (sD):'), disp(sD)
    else
        BMR.applied     = false;
        sD = q; %sD(j) = q;
        rD = p; %rD(j) = p;
        disp('No significant model reduction found.')
    end
else
    % If only one parameter, retain original counts
    sD = q; %sD(j) = q;
    rD = p; %rD(j) = p;
    BMR.applied = false;
    disp('Only one parameter present; no reduction applied.');
end
% Update MDP structure
MDP.d{f}    = sD;
MDP.d0{f}   = rD;
MDP.d_0{f}  = rD;

end

% Score models using Bayesian model reduction
%--------------------------------------------------------------------------
% for i = 1:numel(rD)
%     G    = spm_MDP_log_evidence(qD, pD, pD + rD{i});  % Compute expected reduction in free energy
%     F(i) = sum(G(isfinite(G)));                       % Sum finite elements for each reduced model
% end

% Find any model that has greater evidence than the parent model
%--------------------------------------------------------------------------
% [Fmin, jmin] = min(F);
% if Fmin < -T
%     disp('Fmin < -T')
%     rD_fin = rD{jmin};
% 
%     % Compute data counts and ensure non-negativity
%     data_counts = qD - pD;
%     data_counts(data_counts < 0) = 0;
% 
%     % Compute new posterior concentration parameters under the reduced model
%     sD_fin = rD_fin + data_counts;
% 
%     % Update MDP and BMR structures
%     MDP.d{f}    = sD_fin;
%     MDP.d0{f}   = rD_fin;
%     MDP.d_0{f}  = rD_fin;
% 
%     BMR.yes     = 'yes';
%     BMR.rF      = [BMR.rF {F}];
%     BMR.jmin    = [BMR.jmin jmin];
%     BMR.rFmin   = [BMR.rFmin Fmin];
%     BMR.rD      = [BMR.rD rD_fin];
%     BMR.sD      = [BMR.sD sD_fin];
% 
%     disp('Selected reduced prior (rD_fin):'), disp(rD_fin)
%     disp('Updated posterior (sD_fin):'), disp(sD_fin)
% else
%     disp('No significant model reduction found.')
% end
% end
