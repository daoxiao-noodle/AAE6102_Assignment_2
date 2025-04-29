function skymask()
% GNSS sky occlusion mask visualization (based on skymask_A1_urban.csv)
% The white area represents visible regions, and the gray area represents invisible regions.

% Load occlusion mask data
data = readmatrix('C:\Users\u\Downloads\AAE6102_Assignment2-main\Q2\skymask_A1_urban.csv');
% === GNSS Skyplot Visualization (0° at the top, increasing clockwise) ===

% Step 1: Read the data
%data = readmatrix('skymask_A1_urban.csv');
az_mask = data(:,1);
el_mask = data(:,2);

% Ensure continuity from 0° to 360°
if az_mask(end) ~= 360
    az_mask = [az_mask; 360];
    el_mask = [el_mask; el_mask(1)];
end

% Interpolate to get one value per degree
az_full = 0:359;
el_full = interp1(az_mask, el_mask, az_full, 'linear');

% Step 2: Generate grid and visibility determination
[AZ, EL] = meshgrid(az_full, 0:1:90);
visibility = EL >= repmat(el_full, size(EL,1), 1);

% Polar to Cartesian conversion — ✅ Core correction: clockwise + 0° at the top
theta = deg2rad(-(AZ - 90));  % ✅ Corrected angle mapping
r = 90 - EL;
[X, Y] = pol2cart(theta, r);

% Step 3: Start plotting
figure; 
hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');
xlim([-95 95]); ylim([-95 95]);

% White background circle
rectangle('Position', [-90, -90, 180, 180], 'Curvature', [1 1], 'FaceColor', 'w');

% Occlusion areas (gray patches)
for az = 1:360
    for el = 1:90
        if visibility(el, az) == 0
            X_patch = [X(el,az), X(el,mod(az,360)+1), X(el+1,mod(az,360)+1), X(el+1,az)];
            Y_patch = [Y(el,az), Y(el,mod(az,360)+1), Y(el+1,mod(az,360)+1), Y(el+1,az)];
            patch(X_patch, Y_patch, [0.7 0.7 0.7], 'EdgeColor', 'none');
        end
    end
end

% Step 4: Draw elevation circles (every 15°)
for elevation = 15:15:75
    rectangle('Position', [-elevation, -elevation, 2*elevation, 2*elevation], ...
              'Curvature', [1 1], 'LineStyle', ':', 'EdgeColor', [0.5 0.5 0.5]);
    text(0, elevation, [num2str(90 - elevation), '°'], ...
         'HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w');
end

% Azimuth lines & labels (clockwise)
for az = 0:30:330
    [x_line, y_line] = pol2cart(deg2rad(-(az - 90)), [0, 90]);  % Synchronized clockwise
    plot(x_line, y_line, ':', 'Color', [0.5 0.5 0.5]);

    [x_txt, y_txt] = pol2cart(deg2rad(-(az - 90)), 95);
    switch az
        case 0
            label = '0° (N)';
        case 90
            label = '90° (E)';
        case 180
            label = '180° (S)';
        case 270
            label = '270° (W)';
        otherwise
            label = [num2str(az), '°'];
    end
    text(x_txt, y_txt, label, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

title('skymask');
hold off;
