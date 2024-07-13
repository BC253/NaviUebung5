clear all
close all
clc
addpath("./Funktion")

w_E = [0;0;7.292115e-5];
Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
 %% RTK only
load gnssrtk
rtklat = gnssrtk(:,3);
rtklon = gnssrtk(:,4);
figure
title('rtk')
geoplot(rtklat, rtklon, 'LineWidth', 2, 'Color', 'red')
geobasemap satellite  % 设置为卫星地图
% kmlFileName = 'gnssrtk_path.kml'; 
% kmlwrite(kmlFileName, rtklat, rtklon, 'LineWidth', 2, 'Color', 'red')


%% original GNSS data
opts = detectImportOptions('vn310-gnss.csv'); % 自动检测输入数据的格式，汉字不识别
gnssorginal = readtable('vn310-gnss.csv', opts);%把能识别的部分读出来
% columnNames = gnssorginal.Properties.VariableNames;
% disp(columnNames);  % 显示所有列名
original_lat = table2array(gnssorginal(:,15));
original_lon = table2array(gnssorginal(:,16));
figure
title('original GNSS data')
geoplot(original_lat, original_lon, 'LineWidth', 2, 'Color', 'red')
geobasemap satellite  % 设置为卫星地图
% kmlFileName = 'gnss_original.kml'; 
% kmlwrite(kmlFileName, original_lat, original_lon, 'LineWidth', 2, 'Color', 'red')
%% cleaned GNSS data
opts = detectImportOptions('vn310-gnss.csv'); % 自动检测输入数据的格式，汉字不识别
GNSS.raw = readtable('vn310-gnss.csv', opts);%把能识别的部分读出来
GDOP = table2array(GNSS.raw(:,36));

% 创建逻辑索引，其中GDOP小于2
validIndex = GDOP < 2;

% 使用逻辑索引过滤数据
filtered_lat = original_lat(validIndex);
filtered_lon = original_lon(validIndex);

figure
title('cleaned GNSS data GDOP>2')
geoplot(filtered_lat, filtered_lon, 'LineWidth', 2, 'Color', 'red')
geobasemap satellite  % 设置为卫星地图
% kmlFileName = 'gnss_cleaned.kml'; 
% kmlwrite(kmlFileName, original_lat, original_lon, 'LineWidth', 2, 'Color', 'red')
%% startwert

% 读取IMU数据 ax ay az(m/s^2) gx gy gz(rad/s) 200Hz
optss = detectImportOptions('vn310-imu.csv');
IMU.raw = readtable('vn310-imu.csv', optss);
t_imu = table2array(IMU.raw(:,1));
IMUdata(:,1:3) = table2array(IMU.raw(:,32:34));
IMUdata(:,4:6) = table2array(IMU.raw(:,35:37));
RPY(:,1) = deg2rad(table2array(IMU.raw(:,49)));
RPY(:,2) = deg2rad(table2array(IMU.raw(:,48)));
RPY(:,3) = deg2rad(table2array(IMU.raw(:,47)));
%GNSS 数据
opts = detectImportOptions('vn310-gnss.csv'); % 自动检测输入数据的格式，汉字不识别
GNSS.raw = readtable('vn310-gnss.csv', opts);%把能识别的部分读出来

gnssp_e = table2array(GNSS.raw(:,18:20));  %%初始数据 startwert
gnssp_e_cleaned = table2array(GNSS.raw(validIndex,18:20));
gnssv_e = table2array(GNSS.raw(:,24:26));  %%初始数据 startwert
gnssv_e_cleaned = table2array(GNSS.raw(validIndex,24:26));
gnssv_n = table2array(GNSS.raw(:,21:23));
gnssv_n_cleaned = table2array(GNSS.raw(validIndex,21:23));
gnsslla = table2array(GNSS.raw(:,15:17));
gnsslla = deg2rad(gnsslla(:,1:2));
t_gnss = table2array(GNSS.raw(:,1));
t_gnss_cleaned = table2array(GNSS.raw(validIndex,1));
% GDOP = table2array(GNSS.raw(:,36));
M = length(t_imu);

%DCM t0 berechnen
Cne_0 = C(3,-gnsslla(1,2))*C(2,gnsslla(1,1)+pi/2);
Cbn_0 = C(3,-RPY(1,3))*C(2,-RPY(1,2))*C(1,-RPY(1,1));
Cpe_0 = Cne_0*Cbn_0;
Omega_iep = inv(Cpe_0)*Omega_iee*Cpe_0;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
                                          

%DCM t1 berechnen
Cne_1 = C(3,-gnssp_e(2,2))*C(2,gnssp_e(2,1)+pi/2);
Cbn_1 = C(3,-RPY(2,3))*C(2,-RPY(2,2))*C(1,-RPY(2,1));
Cpe_1 = Cne_1*Cbn_1;
Omega_iep = inv(Cpe_1)*Omega_iee*Cpe_1;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
  
%Quaternion t0 im e-system berechnen 
qt0 = compact(quaternion((Cpe_0),'rotmat','frame'));

%Quaternion t1 im e-system berechnen
qt1 = compact(quaternion((Cpe_1),'rotmat','frame'));
q_e = zeros(M,4);
q_e(1,:) = qt0;
q_e(2,:) = qt1;




%% prepare for KF

updaterate = 1; % 1 oder 5 oder 20, [s]
t_Update = t_gnss(1):updaterate:t_gnss(end);

bias_ap = ones(3,1)*9.81e-5; % initialer Bias Accelerometer
bias_wpip = ones(3,1)*2.42e-6; % initialer Bias Gyro
Rn = diag([1.0 1.0 1.5 0.02 0.02 0.02].^2); % Messrauschen/Unsicherheit Messungen (Pos [m] und Vel [m/s])
Hn = [eye(6), zeros(6,9)]; % Desingmatrix fuer Beobachtungen
xnn = zeros(M, 15); 
xnn(1,:) = ones(15,1)*1e-12; % initialer Zustandsvektor der Fehler [dPos, dVel, dOri, Bias Acc, Bias Gyro]'
Pnn = cell(M, 1);
Pnn{1} = diag([ones(1,3)*1e-1...            
               ones(1,3)*1e-1...           
               deg2rad(ones(1,3)*1e-2)...  
               ones(1,3)*1e-1...           
               deg2rad(ones(1,3)*1)].^2);  % zugehoerige Kovarianzmatrix
xe = zeros(M,3);
xe(1,:) = gnssp_e(1,:);
ve = zeros(M,3);
ve(1,:) = gnssv_e(1,:);
Cpe = Cpe_0;
Rk3_KF = zeros(M,10);
IMUdataKF = zeros(M, 6);
sigma = zeros(M, 15);
%解决时间戳不均匀问题
   tolerance = 1e-5;  
   [isMatch, idx_GNSS] = ismembertol(t_imu, t_gnss, tolerance);
   matchingTimes = t_imu(isMatch);
   matchingIndices = find(isMatch);
for i = 1:length(t_imu)-1    
   
   delta_t = abs(t_imu(i) - t_imu(i+1));


    
   ae = Cpe*IMUdata(i,1:3)';
   Ae = ome2Ome(ae);
   tao = 3000; %更换其他数据时更换
   beta = 1/tao;
   Fa = -beta*eye(3);
   Fw = -beta*eye(3);
   F = [zeros(3,3) eye(3) zeros(3,3) zeros(3,3) zeros(3,3)
     -Omega_iee*Omega_iee -2*Omega_iee Ae Cpe_0 zeros(3)
     zeros(3,3) zeros(3,3) -Omega_iee zeros(3,3) -Cpe_0
     zeros(3,3) zeros(3,3) zeros(3,3) Fa zeros(3,3)
     zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fw];
   varianz_ap = 3.08e-5;    % Varianz Accelerometer Bias Rauschen [m^2/s^4]
   varianz_wpip = 2.424068e-5;  % Varianz Gyro Bias Rauchen [rad^2/s^2]
   Ga = sqrt(varianz_ap*beta)*eye(3); 
   Gw = sqrt(varianz_wpip*beta)*eye(3); 
   G = [zeros(9,6);
     Ga zeros(3,3);
    zeros(3,3) Gw];
   W = 1;
   A = [-F G*W*G'; zeros(size(F)) F']*mean(diff(t_gnss));
   n = length(A)/2;
   B = expm(A);
   Phi = B(n+1:2*n,n+1:2*n)';  % Zustandsuebergangsmatrix
   Q = Phi*B(1:n,n+1:2*n);     % Matrix des Prozessrauschens
      
    if isMatch(i) && idx_GNSS(i) > 0 && (GDOP(idx_GNSS(i))<2) % ja, Beobachtungen vorhanden
        
        zn = [gnssp_e(idx_GNSS(i),:)-xe(i,:) ...  % Positionsdifferenz
              gnssv_e(idx_GNSS(i),:)-ve(i,:)]';   % Geschwindigkeitsdifferenz

        [x, P] = KalmanFilter(xnn(i,:)', Pnn{i} ,Phi, zn, Hn, Q, Rn);
        xnn(i+1,:) = x;
        Pnn{i+1,:} = P;
        sigma(i+1,:) = sqrt(diag(P));

        % geschaetzte Fehler als Korrekturen anbringen (vor Integration
        % (VO06, F7))
        xe(i,:) = xe(i,:) + xnn(i+1,1:3);   % Position
        ve(i,:) = ve(i,:) + xnn(i+1,4:6);   % Geschwindigkeit
            psi = xnn(7:9);
            Psi = [0 -psi(3) psi(2); psi(3) 0 -psi(1); -psi(2) psi(1) 0];
        Cpe = (eye(3) - Psi)*Cpe;     % Orientierung
        q_e(i,:) = compact(quaternion((Cpe),'rotmat','frame'));   % Orientierung
        
        bias_ap = bias_ap+xnn(i+1,10:12)';      % akkumulierter Bias
        bias_wpip = bias_wpip+xnn(i+1,13:end)'; % akkumulierter Bias

        % Fehler zuruecksetzen (Zuschlaege ab jetzt wieder Null, da
        % Korrektur angebracht ist)
        xnn(i+1,1:6) = xnn(1,1:6); 
        
    else
        % Bestimme Fehler fuer den Fall, dass keine Beob. vorliegen
        % (Praediktion Fehlervektor - als Initialwerte fuer nächste Korrektur
        % mit KalmanFilter)
        zn = [];
        [x, P] = KalmanFilter(xnn(i,:)', Pnn{i} ,Phi, zn, Hn, Q, Rn);
        xnn(i+1,:) = x;
        Pnn{i+1,:} = P;
        sigma(i+1,:) = sqrt(diag(P));
    end

    % (korrigierte) Startwerte RungeKutta
    Rk3_KF(i,:) = [xe(i,:),ve(i,:),q_e(i,:)];
    % Anbringen des geschaetzten Bias' an die Messungen
    IMUdataKF(i,1:3) = IMUdata(i,1:3)+bias_ap';
    IMUdataKF(i,4:6) = IMUdata(i,4:6)+bias_wpip';
    if i == 1
        x0 = [IMUdataKF(i, 1:6); IMUdataKF(i, 1:6)];
    else
        x0 = [IMUdataKF(i-1:i, 1:3) IMUdataKF(i-1:i, 4:6)];
    end
    % Integration e-System (RungeKutta 3.Ordnung)
    Rk3_KF(i+1,:) = RungeKutta3(@TimeDerivativePosVelAtt_e_g, Rk3_KF(i,:)', x0', delta_t);

    % Extrahiere Werte aus Integration
    xe(i+1,:) = Rk3_KF(i,1:3);    % Position ECEF [m m m]
    ve(i+1,:) = Rk3_KF(i,4:6);    % Geschwindigkeit ECEF [m/s m/s m/s]
    q_e(i+1,:) = Rk3_KF(i,7:10);         % Quaternion p2e-System
    Cpe  = rotmat(quaternion(q_e(i+1,:)),'frame');% DCM p2e-System
end

x_lla = ecef2lla(xe);

figure
geoplot(filtered_lat, filtered_lon,  'Color', 'red')
hold on
geoplot(x_lla(:,1),x_lla(:,2),'Color','yellow')
legend('gnss_filtered','RK3','REF')
geobasemap satellite



figure
plot3(gnssp_e(:,1),gnssp_e(:,2),gnssp_e(:,3))
hold on
plot3( xe(:,1), xe(:,2), xe(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('3D Plot for all');



figure
plot(gnssp_e(:,2),gnssp_e(:,3))
hold on 
plot(xe(:,2),xe(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('2D Plot for all');



figure
subplot(3, 1, 1);
plot(t_gnss,gnssv_n(:,1));
hold on
plot(t_imu,ve(:,3))
xlabel("[s]");ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in N');

subplot(3, 1, 2);
plot(t_gnss,gnssv_n(:,2));
hold on
plot(t_imu,ve(:,2))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in E');

subplot(3, 1, 3);
plot(t_gnss,gnssv_n(:,3));
hold on
plot(t_imu,ve(:,1))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in D');

%%
figure
plot3(gnssp_e_cleaned(:,1),gnssp_e_cleaned(:,2),gnssp_e_cleaned(:,3))
hold on
plot3( xe(:,1), xe(:,2), xe(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('3D Plot for all');



figure
plot(gnssp_e_cleaned(:,2),gnssp_e_cleaned(:,3))
hold on 
plot(xe(:,2),xe(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('2D Plot for all');



figure
subplot(3, 1, 1);
plot(t_gnss_cleaned,gnssv_n_cleaned(:,1));
hold on
plot(t_imu,ve(:,3))
xlabel("[s]");ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in N');

subplot(3, 1, 2);
plot(t_gnss_cleaned,gnssv_n_cleaned(:,2));
hold on
plot(t_imu,ve(:,2))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in E');

subplot(3, 1, 3);
plot(t_gnss_cleaned,gnssv_n_cleaned(:,3));
hold on
plot(t_imu,ve(:,1))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in D');
