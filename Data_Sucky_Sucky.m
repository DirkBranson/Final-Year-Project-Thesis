
clear all
clc
close all
    if ~isempty(instrfind)
     fclose(instrfind);
      delete(instrfind);
    end
s = serial('COM12', 'BaudRate',115200);
 fopen(s);
 
Table_Header = {'Accel_X','Accel_Y','Accel_Z','Gyro_X','Gyro_Y','Gyro_Z','MagX_uT','MagY_uT','MagZ_uT','Temp_C','Time'};
count = 0;
count_final = 500;
Accel_X = zeros(1,count_final);
Accel_Y = zeros(1,count_final);
Accel_Z = zeros(1,count_final);
Gyro_X = zeros(1,count_final);
Gyro_Y = zeros(1,count_final);
Gyro_Z = zeros(1,count_final);
MagX_uT = zeros(1,count_final);
MagY_uT = zeros(1,count_final);
MagZ_uT = zeros(1,count_final);

Data_Initial = split(fscanf(s),',');
while count <= count_final
  count = count + 1
        Data{count} = split(fscanf(s),',');
        Accel_X(1,count) = str2num(cell2mat(Data{count}(2,1))) - str2num(cell2mat(Data{1}(2,1)))
        Accel_Y(count) = str2num(cell2mat(Data{count}(3,1))) - str2num(cell2mat(Data{1}(3,1)))
        Accel_Z(count) = str2num(cell2mat(Data{count}(4,1))) - str2num(cell2mat(Data{1}(4,1)))
        Gyro_X(count) = str2num(cell2mat(Data{count}(5,1))) - str2num(cell2mat(Data{1}(5,1)));
        Gyro_Y(count) = str2num(cell2mat(Data{count}(6,1))) - str2num(cell2mat(Data{1}(6,1)));
        Gyro_Z(count) = str2num(cell2mat(Data{count}(7,1))) - str2num(cell2mat(Data{1}(7,1)));
        MagX_uT(count) = str2num(cell2mat(Data{count}(8,1))) - str2num(cell2mat(Data{1}(8,1)));
        MagY_uT(count) = str2num(cell2mat(Data{count}(9,1))) - str2num(cell2mat(Data{1}(9,1)));
        MagZ_uT(count) = str2num(cell2mat(Data{count}(10,1))) - str2num(cell2mat(Data{1}(10,1)));
        
        %%Countinuous timer modification
              if count == 1 
            Time(count) = str2num(cell2mat(Data{count}(11,1)));
        else 
                Time(count) = Time(count - 1) + str2num(cell2mat(Data{count}(11,1)));
              end
        
        Accel_Mag(count) = (Accel_X(count)^2 + Accel_Y(count)^2 + Accel_Z(count)^2)^0.5;
        Gyro_Mag(count) = (Gyro_X(count)^2 + Gyro_Y(count)^2 + Gyro_Z(count)^2)^0.5;
        Mag_Mag(count) = (MagX_uT(count)^2 + MagY_uT(count)^2 + MagZ_uT(count)^2)^0.5;
end
fclose(s);

%Map these to MagMag

figure
colormap jet
xlabel({'Gyro_X'})
ylabel({'Gyro_Y'})
zlabel({'Gyro_Z'})
scatter3(Gyro_X,Gyro_Y,Gyro_Z,[],Mag_Mag)
caxis([0 10000])
colorbar

figure
colormap jet
xlabel({'Gyro_X'})
ylabel({'Gyro_Y'})
zlabel({'Time/s'})
scatter3(Gyro_X,Gyro_Y,Time,[],Mag_Mag)
caxis([0 10000])
colorbar

figure
colormap jet
scatter3(Accel_X,Accel_Y,Time,[],Mag_Mag,'x')
xlabel({'Accel_X'})
ylabel({'Accel_Y'})
zlabel({'Time/s'})
xlim([-10 10])
ylim([-10 10])
zlim([0 1.5])
caxis([0 10000])
colorbar

[Accel_X_Noise_Peaks,Accel_X_Noise_Peaks_Locs] = findpeaks(Accel_X(1:50),'SortStr','descend','NPeaks',5)
[Accel_Y_Noise_Peaks,Accel_Y_Noise_Peaks_Locs] = findpeaks(Accel_Y(1:50),'SortStr','descend','NPeaks',5)
[Accel_Z_Noise_Peaks,Accel_Z_Noise_Peaks_Locs] = findpeaks(Accel_Z(1:50),'SortStr','descend','NPeaks',5)
%Moving average filter
windowSize = 20; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Accel_X_Filter = filter(b,a,Accel_X);
Accel_Y_Filter = filter(b,a,Accel_Y);
Accel_Z_Filter = filter(b,a,Accel_Z);

for j = 1:count
    if abs(Accel_X_Filter(j)) <= 0.1
    Accel_X_Filter2(j) = 0;
    else
        Accel_X_Filter2(j) = Accel_X(j)
    end
end
for j = 1:count
    if abs(Accel_Y_Filter(j)) <= 0.1
    else
        Accel_Y_Filter2(j) = Accel_Y(j)
    end
end
for j = 1:count
    if abs(Accel_Z_Filter(j)) <= 0.1
    Accel_Z_Filter2(j) = 0;
    else
        Accel_Z_Filter2(j) = Accel_Z(j)
    end
end

Accel_X0 = 0;
Velo_X0 = 0;

Accel_wrt_Time = interp1(Time,Accel_X_Filter2,'linear','pp');
Accel_func = @(seconds) ppval(seconds,Accel_wrt_Time);
Accel_func(4)
Velo_func = @(seconds) integral(Accel_func,Time(1),seconds) + Accel_X0;
for i = 1:length(Time)
    Velo_X(i) = Velo_func(i);
end

Velo_wrt_Time = interp1(Time,Velo_X,'linear','pp');
Velo_func = @(seconds) ppval(seconds,Velo_wrt_Time);
Velo_func(4)
Disp_func = @(seconds) integral(Velo_func,Time(1),seconds) + Velo_X0;
for i = 1:length(Time)
    Disp_X(i) = Disp_func(i);
end
Disp_wrt_Time = interp1(Time,Disp_X,'linear','pp');


figure
hold on
plot(Time,Disp_X,'g')
plot(Time,Velo_X,'r')
plot(Time,Accel_X,'b')
plot(Time,Accel_X_Filter,'k')
plot(Time,Accel_X_Filter2,'c')
 ylim([ -10 10])
 xlim([0 260])

figure
plot(Time,(Mag_Mag/1000000))
title('Magnetic Flux Density over Time in Induction Motor')
ylabel('Magnetic Flux Density')
xlabel('Time')




