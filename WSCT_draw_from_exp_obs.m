function [MDP] = WSCT_draw_from_exp_obs(MDP, OBS)
for i=1:length(MDP)
    for feature=1:3
        feature_vec = zeros(numel(MDP(1).D{1}), 1);
        feature_vec((OBS{i}(feature, 1))) = 1;
        MDP(i).D{feature} = feature_vec;
    end
end
end