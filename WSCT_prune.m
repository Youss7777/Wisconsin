function [sD_fin,rD_fin] = WSCT_prune(qD,pD,x)
% FORMDT [sD,rD] = WSCT_prune(qD,pD,f,x,T,m)
% qD - posterior expectations
% pD - prior expectations
% f  - hidden factor to integrate over [defult: 0]
% x  - prior counts [default: 8]
% T  - threshold for Bayesian model reduction [default: three]
%
% sD - reduced posterior expectations
% rD - reduced prior expectations
%__________________________________________________________________________



% model space: additional concentration parameters (i.e., precision)
%--------------------------------------------------------------------------

n = numel(pD);
motifs = ff2n(n);           % all combinations of binary vectors of dim n
motifs(1, : ) = [];         % remove all zeros vectors
motifs(end, : ) = [];       % remove all ones vectors

for mot=1:size(motifs, 1)
    rD{mot} = motifs(mot, :)'*x;
end


% score models using Bayesian model reduction
%--------------------------------------------------------------------------
for i = 1:numel(rD)
    G    = spm_MDP_log_evidence(qD,pD,pD + rD{i});  % compute expected reduction in free energy
    F(i) = sum(G(isfinite(G)));                     % for each reduced model
end

% find any model that has greater evidence than the parent model
%--------------------------------------------------------------------------
[Fmin,jmin] = min(F);
if Fmin < 0
    rD_fin = rD{jmin};
else
    rD_fin = spm_zeros(pD);
end

sD_fin = qD + rD_fin - pD;


return

