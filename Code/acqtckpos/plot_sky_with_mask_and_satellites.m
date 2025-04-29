function plot_sky_with_mask_and_satellites(navSolutionsWLS, sv, skymask_csv_path)

% === Step 1: 读取遮挡掩膜数据 ===
data = readmatrix(skymask_csv_path);
az_mask = data(:,1);
el_mask = data(:,2);

if az_mask(end) ~= 360
    az_mask = [az_mask; 360];
    el_mask = [el_mask; el_mask(1)];
end

az_full = 0:359;
el_full = interp1(az_mask, el_mask, az_full, 'linear');

% 网格生成
[AZ, EL] = meshgrid(az_full, 0:1:90);
visibility = EL >= repmat(el_full, size(EL,1), 1);

% 极坐标转直角坐标
theta = deg2rad(-(AZ - 90));  
r = 90 - EL;
[X, Y] = pol2cart(theta, r);

% === Step 2: 开始绘图 ===
figure; hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');
xlim([-95 95]); ylim([-95 95]);

% 画白色背景圆
rectangle('Position',[-90,-90,180,180],'Curvature',[1,1],'FaceColor','w');

% 画遮挡区域（灰色）
for az = 1:360
    for el = 1:90
        if visibility(el, az) == 0
            X_patch = [X(el,az), X(el,mod(az,360)+1), X(el+1,mod(az,360)+1), X(el+1,az)];
            Y_patch = [Y(el,az), Y(el,mod(az,360)+1), Y(el+1,mod(az,360)+1), Y(el+1,az)];
            patch(X_patch, Y_patch, [0.7 0.7 0.7], 'EdgeColor', 'none');
        end
    end
end

% 绘制仰角圈（每15°）
for elevation = 15:15:75
    rectangle('Position',[-elevation,-elevation,2*elevation,2*elevation],...
        'Curvature',[1,1],'LineStyle',':','EdgeColor',[0.5 0.5 0.5]);
    text(0,elevation,[num2str(90-elevation),'°'],'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','BackgroundColor','w');
end

% 绘制方位角线（每30°）
for az = 0:30:330
    [x_line, y_line] = pol2cart(deg2rad(-(az - 90)), [0,90]);
    plot(x_line, y_line, ':', 'Color', [0.5 0.5 0.5]);
    
    [x_txt, y_txt] = pol2cart(deg2rad(-(az - 90)), 95);
    switch az
        case 0, label = '0° (N)';
        case 90, label = '90° (E)';
        case 180, label = '180° (S)';
        case 270, label = '270° (W)';
        otherwise, label = [num2str(az), '°'];
    end
    text(x_txt, y_txt, label, 'HorizontalAlignment','center','VerticalAlignment','middle');
end

% === Step 3: 画卫星点 ===
az_sats = navSolutionsWLS.satAZ(1, :);  % 拿最后一个历元，所有卫星
el_sats = navSolutionsWLS.satEA(1, :);
% az_sats = navSolutionsWLS.satAZ(end, :);  % 拿最后一个历元，所有卫星
% el_sats = navSolutionsWLS.satEA(end, :);
prn_sats = sv;                            % PRN号就是sv

% az_sats 和 el_sats 的顺序已经和 sv一致了


for i = 1:length(prn_sats)
    % 计算极坐标
    r_sat = 90 - el_sats(i);
    theta_sat = deg2rad(-(az_sats(i) - 90));
    [x_sat, y_sat] = pol2cart(theta_sat, r_sat);

    % 判断是否可见（用skymask数据）
    az_idx = mod(round(az_sats(i)),360) + 1; 
    el_idx = min(90, max(0, round(el_sats(i))));
    visible = el_sats(i) >= el_full(az_idx);

    % 根据可见性设定颜色
    if visible
        color = 'g'; % 绿色可见
    else
        color = 'r'; % 红色遮挡
    end

    % 画卫星点
    plot(x_sat, y_sat, 'o', 'MarkerSize',8, 'MarkerFaceColor', color, 'MarkerEdgeColor', color);

    % 标注PRN号
    text(x_sat, y_sat, sprintf('%d', prn_sats(i)), ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'Color','black', 'FontSize',8, 'FontWeight','bold');
end

title('Skyplot with Skymask and Satellites');
hold off;
