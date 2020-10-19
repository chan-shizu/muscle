clc;
clear;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
path_tracking_arm_surface = ["D:\Documents\作業まとめ\2020_0923\result.csv"];
tracking_arm_surface = readmatrix(path_tracking_arm_surface);
tracking_arm_surface = tracking_arm_surface(10:end,:)

tracking_arm_surface_6 = tracking_arm_surface(:,6)
tracking_arm_surface_7 = tracking_arm_surface(:,7)
tracking_arm_surface_14 = tracking_arm_surface(:,14)
tracking_arm_surface_15 = tracking_arm_surface(:,15)
tracking_arm_surface(:,6) = tracking_arm_surface_7
tracking_arm_surface(:,7) = tracking_arm_surface_6
tracking_arm_surface(:,14) = tracking_arm_surface_15
tracking_arm_surface(:,15) = tracking_arm_surface_14

tracking_arm_surface = tracking_arm_surface(10:end,:)
size_tracking = size(tracking_arm_surface);
tracking_arm_surface_x = tracking_arm_surface(:,2:7);
tracking_arm_surface_y = 900 - tracking_arm_surface(:,11:16);

for i=1:size_tracking(1)
    plot(tracking_arm_surface_x(i,:),tracking_arm_surface_y(i,:));
    axis([50 650 500 900])
     Frame(i) = getframe(1);
end

v = VideoWriter('speedup');
videoName = ["D:\pictures\result_graph"];
v = VideoWriter(videoName);
v.Quality=50;
v.FrameRate = 10

open(v);
writeVideo(v,Frame);
close(v);
