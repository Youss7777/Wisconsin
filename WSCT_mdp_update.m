function MDP = spm_MDP_update(MDP,OUT)
% FORMAT MDP = spm_MDP_update(MDP,OUT)
% moves Dirichlet parameters from OUT to MDP
% MDP - structure array (new)
% OUT - structure array (old)
%__________________________________________________________________________

% check for concentration parameters at this level
%--------------------------------------------------------------------------
try,  MDP.a = OUT.a; end
try,  MDP.b = OUT.b; end
try,  MDP.c = OUT.c; end
try,  MDP.d = OUT.d; end
try,  MDP.d_0 = OUT.d_0; end        % added this for BMR
try,  MDP.d0 = OUT.d0; end          % and this
try,  MDP.e = OUT.e; end

% check for concentration parameters at nested levels
%--------------------------------------------------------------------------
try,  MDP.MDP(1).a = OUT.mdp(end).a; end
try,  MDP.MDP(1).b = OUT.mdp(end).b; end
try,  MDP.MDP(1).c = OUT.mdp(end).c; end
try,  MDP.MDP(1).d = OUT.mdp(end).d; end
try,  MDP.MDP(1).e = OUT.mdp(end).e; end

return


function L = spm_MDP_VB_VOX(MDP,L,t)
% FORMAT L = spm_MDP_VB_VOX(MDP,L,t)
% returns likelihoods from voice recognition (and articulates responses)
% MDP - structure array
% L   - predictive prior over outcomes
% t   - current trial
%
% L   - likelihood of lexical and prosody outcomes
%
% this subroutine determines who is currently generating auditory output
% and produces synthetic speech - or uses the current audio recorder object
% to evaluate the likelihood of the next word
%__________________________________________________________________________


% check for VOX structure
%--------------------------------------------------------------------------
global VOX
global TRAIN
if ~isstruct(VOX), load VOX; VOX.RAND = 0; end
if isempty(TRAIN), TRAIN = 0;              end
if t == 1,         pause(1);               end


if ~isfield(VOX,'msg')
    
    % prepare useful fields in (global) VOX structure
    %----------------------------------------------------------------------
    Data    = imread('recording','png');
    VOX.msg = msgbox('Recording','','custom',Data);
    set(VOX.msg,'Visible','off'), drawnow
    
    % indices of words in lexicon and inferred prosody states
    %----------------------------------------------------------------------
    VOX.io  = spm_voice_i(MDP.label.outcome{1});
    VOX.ip  = find(ismember({VOX.PRO.str},{'amp','dur','Tf','p0','p1','p2'}));
    
    % check for audio recorder
    %----------------------------------------------------------------------
    if ~isfield(VOX,'audio')
        VOX.audio  = audiorecorder(22050,16,1);
    end
    
end

if MDP.VOX == 0 || MDP.VOX == 1
    
    % Agent: computer
    %----------------------------------------------------------------------
    str = MDP.label.outcome{1}(MDP.o(1,1:t));
    fprintf('%i: %s\n',MDP.VOX, str{t});
    i   = ismember(str,' ');
    str = str(~i);
    eof = sum(i);
    
    % if this is the end of a sentence
    %----------------------------------------------------------------------
    if eof == 1 || (~eof && t == MDP.T)
        
        % get lexical and prosody
        %------------------------------------------------------------------
        lexical = spm_voice_i(str);
        prosody = VOX.prosody(:,lexical);
        if MDP.VOX == 0
            speaker = [12;12];
        else
            speaker = [3; 3 ];
        end
        
        % add prosody and articulate
        %------------------------------------------------------------------
        prosody(VOX.ip,:) = MDP.o(2:end,1:numel(str));
        prosody(1,:)      = min(8,prosody(1,:) + 2);
        spm_voice_speak(lexical,prosody,speaker);
        
        
        % TRAIN: prompt for prosody
        %------------------------------------------------------------------
        if TRAIN
            VOX.mute  = 0;
            VOX.depth = 1;
            [i,P]     = spm_voice_i(str);
            prosody   = [];
            while size(prosody,2) ~= numel(str)
                clc
                disp('Please repeat:'), disp(str)
                [SEG,W,prosody] = spm_voice_read(VOX.audio,P);
            end
            clc, disp('Thank you')
            
            % prompt for prosody
            %--------------------------------------------------------------
            for i = 1:numel(VOX.ip)
                for j = 1:size(P,2)
                    L{i + 1,j} = sparse(prosody(VOX.ip(i),j),1,1,8,1);
                end
            end
            
            % uniform priors for spce  (' ')
            %--------------------------------------------------------------
            L{i + 1,t} = ones(8,1)/8;
            
        end
        
        
    end

    
elseif MDP.VOX == 2
    
    % user
    %----------------------------------------------------------------------
    if t == 1
        VOX.IT = 1;
        stop(VOX.audio)
        record(VOX.audio,8);
        set(VOX.msg,'Visible','on')
        pause(1);
        set(VOX.msg,'Visible','off')
        
        % toggle to see spectral envelope
        %----------------------------------------------------------------------
        VOX.onsets = 0;
        
    end
    
    % get prior over outcomes and synchronise with Lexicon
    %----------------------------------------------------------------------
    io  = VOX.io;                            % indices words in lexicon
    ip  = VOX.ip;                            % indices of prosody
    no  = numel(io);                         % number of outcomes
    nw  = numel(VOX.LEX);                    % number of words in lexicon
    nk  = size(L,2) - t + 1;                 % number of predictions
    P   = zeros(nw,nk);                      % prior over lexicon
    for k = 1:nk
        for i = 1:no
            j = io(i);
            if j
                P(j,k) = L{1,t + k - 1}(i);
            end
        end
    end
      
    % deep segmentation: check for last word (P(:,2) = 0)
    %----------------------------------------------------------------------
    VOX.LL = -128;
    VOX.LW = 0;
    if size(P,2) > 1
        if any(P(:,2))
            if any(P(:,2) < (1 - 1/8))
                VOX.LL = 4;                  % there may be no next word
                VOX.LW = 0;
            else
                VOX.LL = -128;               % there is a subsequent word
                VOX.LW = 0;                  % and this is not the last
            end
        else
            P      = P(:,1);
            VOX.LW = 1;                      % this is the last word
        end
    end
    
    % or direct segmentation (comment out to suppress)
    %----------------------------------------------------------------------
    P  = P(:,1);
    
    % get likelihood of discernible words
    %----------------------------------------------------------------------
    P      = P > 1/128;
    if any(P(:,1))
        
        % log likelihoods
        %------------------------------------------------------------------
        O  = spm_voice_get_word(VOX.audio,bsxfun(@rdivide,P,sum(P)));
        
        % check for end of sentence 
        %------------------------------------------------------------------
        try
            A = spm_softmax(O{2}(:,1));    % P(lowest amplitude)
        catch
            A = 1;
        end
        if isempty(O) || A(1) > 1/2
            
            % end of sentence or indiscernible word
            %--------------------------------------------------------------
            L{1,t}( ~io) = 1;
            L{1,t}(~~io) = 0;
            L{1,t}       = L{1,t}/sum(L{1,t});
            
            % prosody likelihoods
            %--------------------------------------------------------------
            for g = 2:numel(L)
                L{g,t} = spm_softmax(spm_zeros(L{g,t}));
            end
            
        else
            
            % lexical likelihoods
            %--------------------------------------------------------------
            LL    = zeros(no,1);
            for i = 1:no
                j = io(i);
                if j
                    LL(i) = O{1}(j);
                end
            end
            L{1,t}  = spm_softmax(LL);
                        
            % prosody likelihoods
            %--------------------------------------------------------------
            for g = 2:numel(L)
                L{g,t} = spm_softmax(O{2}(:,ip(g - 1)));
            end

        end

    end % discernible words
    
    % display word
    %----------------------------------------------------------------------
    [d,w]  = max(L{1,t});
    fprintf('%i: %s\n',MDP.VOX, MDP.label.outcome{1}{w});
        

end