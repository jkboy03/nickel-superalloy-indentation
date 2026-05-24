%% IPDF plots of indentation modulus (E_ind) and yield (Y_ind)
%
% Generates inverse pole distribution function (IPDF) plots colored by the
% measured E_ind and Y_ind, with the z-axis as the reference direction.
% Requires MTEX (https://mtex-toolbox.github.io) on the MATLAB path.
%
% Cleaned-up version of the original plot_ipdf.m. Single CSV read, no dead
% code, relative paths, and configurable color ranges and save options.

%% Startup
startup_mtex                    % only needed once per MATLAB session
clear all
close all
clc

%% ----- Configuration ---------------------------------------------------
% EDIT THIS if you move the project: absolute path to the Bayesian Paper
% folder. Hardcoded as an absolute path so the script works no matter how
% you launch it (MATLAB editor "Run Section" sometimes copies the file to
% a temp dir, breaking mfilename('fullpath')).
project_root = '';
data_file    = fullfile(project_root, 'data', 'experimental', 'ar_final.csv');

% Sanity check
if ~exist(data_file, 'file')
    error('plot_ipdf:DataFileNotFound', ...
          ['Could not find data file at:\n  %s\n\nEdit ' ...
           '`project_root` near the top of this script to match where ' ...
           'your Bayesian Paper folder is on disk.'], data_file);
end

% Color scale ranges. Set to [] to let MATLAB auto-range.
clim_E_ind  = [242, 267];       % GPa (indentation modulus)
clim_Y_ind  = [2.0,  2.47];     % GPa (yield / hardness)

marker_size = 20;
font_size   = 20;

% Set true to write PNGs into <project_root>/figures/ alongside the Python output.
save_figs   = false;
fig_dir     = fullfile(project_root, 'figures');

%% ----- Load experimental data ------------------------------------------
T = readtable(data_file);
fprintf('Loaded %d indentation tests from %s\n', height(T), data_file);

% Crystal symmetry and per-test orientations from Euler angles (radians)
cs  = crystalSymmetry('m-3m');                                  % FCC
ori = orientation.byEuler(T.phi1, T.PHI, T.phi2, cs);

%% ----- IPF color key (reference plot) ----------------------------------
figure;
plot(ipfColorKey(cs));
title('IPF color key (m-3m, FCC)');

%% ----- IPDF colored by E_ind -------------------------------------------
figure;
for i = 1:height(T)
    plotIPDF(ori(i), T.Eind(i), zvector, 'antipodal', ...
             'MarkerEdgeColor', [0 0 0], ...
             'MarkerSize', marker_size);
    hold on
end
if ~isempty(clim_E_ind),  caxis(clim_E_ind);  end
colormap(coolwarm(256));
cb = colorbar;
cb.Label.String = 'E_{ind} (GPa)';
ax = gca;  ax.FontSize = font_size;  ax.Color = 'w';

if save_figs
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    print(fullfile(fig_dir, 'ipdf_Eind.png'), '-dpng', '-r300');
end

%% ----- IPDF colored by Y_ind -------------------------------------------
figure;
for i = 1:height(T)
    plotIPDF(ori(i), T.Yind(i), zvector, 'antipodal', ...
             'MarkerEdgeColor', [0 0 0], ...
             'MarkerSize', marker_size);
    hold on
end
if ~isempty(clim_Y_ind),  caxis(clim_Y_ind);  end
colormap(coolwarm(256));
cb = colorbar;
cb.Label.String = 'Y_{ind} (GPa)';
ax = gca;  ax.FontSize = font_size;  ax.Color = 'w';

if save_figs
    print(fullfile(fig_dir, 'ipdf_Yind.png'), '-dpng', '-r300');
end
