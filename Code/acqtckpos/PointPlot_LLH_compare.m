function PointPlot_LLH_compare(navSolutions1,navSolutions2, navSolutions3, navSolutions4, navSolutions5, iniPos)
%function PointPlot_LLH_compare(navSolutions1,navSolutions2, iniPos)
%PointPlot_LLH_compare(navSolutionsWLS,navSolutionsWLS2,navSolutionsWLS4,navSolutionsWLS_TASK5, solu.iniPos)
% Function to plot GPS estimated positions on a geographic map
% Inputs:
%   navSolutions - Structure containing estimated GPS positions
%   iniPos - Reference position [latitude, longitude, altitude] in radians

% Extract latitude, longitude, and altitude from estimated positions
latitudes = navSolutions1.usrPosLLH(:,1); % Estimated latitude
longitudes = navSolutions1.usrPosLLH(:,2); % Estimated longitude
altitudes = navSolutions1.usrPosLLH(:,3); % Estimated altitude

latitudes2 = navSolutions2.usrPosLLH(:,1); % Estimated latitude
longitudes2 = navSolutions2.usrPosLLH(:,2); % Estimated longitude
altitudes2 = navSolutions2.usrPosLLH(:,3); % Estimated altitude

latitudes3 = navSolutions3.usrPosLLH(:,1); % Estimated latitude
longitudes3 = navSolutions3.usrPosLLH(:,2); % Estimated longitude
altitudes3 = navSolutions3.usrPosLLH(:,3); % Estimated altitude

latitudes4 = navSolutions4.usrPosLLH(:,1); % Estimated latitude
longitudes4 = navSolutions4.usrPosLLH(:,2); % Estimated longitude
altitudes4 = navSolutions4.usrPosLLH(:,3); % Estimated altitude

latitudes5 = navSolutions5.usrPosLLH(:,1); % Estimated latitude
longitudes5 = navSolutions5.usrPosLLH(:,2); % Estimated longitude
altitudes5 = navSolutions5.usrPosLLH(:,3); % Estimated altitude


% Convert reference position from radians to degrees
RefLat = iniPos(1) * 180 / pi;
RefLon = iniPos(2) * 180 / pi;
RefAlt = iniPos(3) * 180 / pi;

% Create geographic map
figure;
geoscatter(latitudes, longitudes, 20, 'blue', 'filled', 'DisplayName', 'Estimated Position');
hold on;
geoscatter(latitudes2, longitudes2, 20, 'g', 'filled', 'DisplayName', 'Estimated Position');
% Plot reference position as a red dot
hold on;
geoscatter(latitudes3, longitudes3, 20, 'yellow', 'filled', 'DisplayName', 'Estimated Position');
hold on;
geoscatter(latitudes4, longitudes4, 20, 'cyan', 'filled', 'DisplayName', 'Estimated Position');
hold on;
geoscatter(latitudes5, longitudes5, 20, 'm', 'filled', 'DisplayName', 'Estimated Position');
hold on;

geoscatter(RefLat, RefLon, 200, 'r', 'filled', 'DisplayName', 'Ground Truth');

% Set satellite basemap
geobasemap('satellite'); % Options: 'satellite', 'topographic', 'streets'

% Set title
title('Position Results on Geographic Map');

% Display legend
%legend('WLS','WLS+RAIM', 'best');
legend('weight=1','weight=0.8','weight=0.5','weight=0.3', 'weight=0','best');

hold off;

% ==== 计算每条轨迹的RMSE ====

% 将真值位置iniPos转成[经度, 纬度, 高度]，单位是degrees（要和usrPosLLH一致）
%truth = [iniPos(1)*180/pi, iniPos(2)*180/pi, iniPos(3)]; 
truth = [0, 0, 0]; 
% 每一条轨迹的position matrix
pos1 = [navSolutions1.usrPosENU(:,1), navSolutions1.usrPosENU(:,2), navSolutions1.usrPosENU(:,3)];
pos2 = [navSolutions2.usrPosENU(:,1), navSolutions2.usrPosENU(:,2), navSolutions2.usrPosENU(:,3)];
pos3 = [navSolutions3.usrPosENU(:,1), navSolutions3.usrPosENU(:,2), navSolutions3.usrPosENU(:,3)];
pos4 = [navSolutions4.usrPosENU(:,1), navSolutions4.usrPosENU(:,2), navSolutions4.usrPosENU(:,3)];
pos5 = [navSolutions5.usrPosENU(:,1), navSolutions5.usrPosENU(:,2), navSolutions5.usrPosENU(:,3)];

% 计算每条轨迹到真值的位置误差（欧式距离），再取RMSE
rmse1 = sqrt(mean(sum((pos1 - truth).^2, 2)));
rmse2 = sqrt(mean(sum((pos2 - truth).^2, 2)));
rmse3 = sqrt(mean(sum((pos3 - truth).^2, 2)));
rmse4 = sqrt(mean(sum((pos4 - truth).^2, 2)));
rmse5 = sqrt(mean(sum((pos5 - truth).^2, 2)));

% 输出结果
fprintf('RMSE of navSolutions1: %.4f meters\n', rmse1);
fprintf('RMSE of navSolutions2: %.4f meters\n', rmse2);
fprintf('RMSE of navSolutions3: %.4f meters\n', rmse3);
fprintf('RMSE of navSolutions4: %.4f meters\n', rmse4);
fprintf('RMSE of navSolutions5: %.4f meters\n', rmse5);

end
