function [F,Fu,Fs,Fq,Fg,Fd, p_best_policy, prec] = WSCT_spm_MDP_F(MDP, f_state)
% auxiliary function for retrieving free energy and its components
% FORMAT [F,Fu,Fs,Fq,Fg,Fa] = spm_MDP_F(MDP)
%
% F   - total free energy
% Fu  - confidence
% Fs  - free energy of states
% Fq  - free energy of policies
% Fg  - free energy of precision
% Fd  - free energy of parameters
%
% If MDP is a cell array, the free actions are turned (summed over time),
% otherwise, the free energies are turned over time
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_F.m 6811 2016-06-17 09:55:47Z karl $


% evaluate free action
%==========================================================================
m = numel(MDP);
if m > 1
    for i = 1:numel(MDP)
        
        % free action due to states and policies
        %------------------------------------------------------------------
        [f,fu,fs,fq,fg] = WSCT_spm_MDP_F(MDP(i));
        F(i) = sum(f);
        Fu(i) = sum(fu);
        Fs(i) = sum(fs);
        Fq(i) = sum(fq);
        Fg(i) = sum(fg);
        Fd(i) = MDP(i).Fd(f_state);
        act_prob(:,i) = MDP(i).P(:, :, :, :, :, :, 2);
        p_best_policy(i) = max(act_prob(:,i));
        prec(:,i) = MDP(i).w(3);

        % free energy due to parameters
        %------------------------------------------------------------------
        try
            Fd = spm_cat({MDP.Fd});
        catch
            Fd = 0;
        end
    end
    return
    
else
    
    % evaluate free energy
    %======================================================================
    pg  = 1;                                    % prior precision
    qg  = MDP.w;                                % posterior precision
    pu  = spm_softmax(MDP.G*diag(qg));          % prior policies
    qu  = spm_softmax(MDP.F + MDP.G*diag(qg));  % posterior policies
    
    Fu  =  sum(qu.*log(qu));                    % confidence
    Fs  = -sum(qu.*MDP.F);                      % free energy of states
    Fq  = -sum(qu.*log(pu));                    % free energy of policies
    Fg  = qg/pg - log(qg);                      % free energy of precision
    Fd  = [];
    F   = Fs + Fu + Fq + Fg;                    % total free energy
    act_prob = MDP.P(:, :, :, :, :, :, 2);
    p_best_policy = max(act_prob);
    
end