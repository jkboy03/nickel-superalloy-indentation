function cmap = coolwarm(n)
%COOLWARM  matplotlib's 'coolwarm' diverging colormap for MATLAB.
%
%   CMAP = COOLWARM(N) returns an Nx3 RGB colormap matching matplotlib's
%   coolwarm — a smooth blue/gray/red diverging map (Moreland, 2009).
%   Default N is 256.
%
%   Use it like any other colormap:
%       colormap(coolwarm(256));
%       colormap(coolwarm);          % default 256 levels
%
%   See also: jet, parula, viridis.

    if nargin < 1 || isempty(n)
        n = 256;
    end

    % Anchor colors copied directly from matplotlib's _coolwarm_data.
    % 33 stops evenly spaced from cool (deep blue) to warm (deep red).
    anchors = [
        0.2298057,  0.298717966, 0.753683153;
        0.26623388, 0.353094838, 0.801466763;
        0.30386891, 0.406535296, 0.84495867;
        0.342804478,0.458757618, 0.883725899;
        0.38301334, 0.50941904,  0.917387822;
        0.424369608,0.558148092, 0.945619588;
        0.46666708, 0.604562568, 0.968154911;
        0.509635204,0.648280772, 0.98478814;
        0.552953156,0.688929332, 0.995375608;
        0.596262162,0.726149107, 0.999836203;
        0.639176211,0.759599947, 0.998151185;
        0.681291281,0.788964712, 0.990363227;
        0.722193294,0.813952739, 0.976574709;
        0.761464949,0.834302879, 0.956945269;
        0.798691636,0.849786142, 0.931688648;
        0.833466556,0.860207984, 0.901068838;
        0.865395197,0.86541021,  0.865395561;
        0.897787179,0.848937047, 0.820880546;
        0.924127593,0.827384882, 0.774508472;
        0.944468518,0.800927443, 0.726736146;
        0.958852946,0.769767752, 0.678007945;
        0.96732803, 0.734132809, 0.628751763;
        0.969954137,0.694266682, 0.579375448;
        0.966811177,0.650421156, 0.530263762;
        0.958003065,0.602842431, 0.481775914;
        0.943660866,0.551750968, 0.434243684;
        0.923944917,0.49730856,  0.387970225;
        0.89904617, 0.439559467, 0.343229596;
        0.869186849,0.378313092, 0.300267182;
        0.834620542,0.312874446, 0.259301199;
        0.795631745,0.24128379,  0.220525627;
        0.752534934,0.157246067, 0.184115123;
        0.705673158,0.01555616,  0.150232812;
    ];

    % Smooth interpolation (pchip avoids overshoot at the anchor points).
    x_anchors = linspace(0, 1, size(anchors, 1));
    x_out     = linspace(0, 1, n);
    cmap = interp1(x_anchors, anchors, x_out, 'pchip');
    cmap = max(0, min(1, cmap));     % clip to valid [0, 1] range
end
