%プログラムの高速化のために事前割り当て
datai_x = zeros(timeNum(1),y*t*h);
datai_y = datai_x;
datai_z = datai_x;
con_x = datai_x;
con_y = datai_x;
con_z = datai_x;
Fm_xFm_y = datai_x;
Fm_z = datai_x;
Fn_x = datai_x;
Fn_y = datai_x;
Fn_z = datai_x;
Fsl_x = datai_x;
Fsl_y = datai_x;
Fsl_z = datai_x;
Fsr_x = datai_x;
Fsr_y = datai_x;
Fsr_z = datai_x;
FsSum_x = datai_x;
FsSum_y = datai_x;
FsSum_z = datai_x;
Fv_x = datai_x;
Fv_y = datai_x;
Fv_z = datai_x;
vis_x = datai_x;
vis_y = datai_x;
vis_z = datai_x;
vn_x = datai_x;
vn_y = datai_x;
vn_z = datai_x;
Mr_x = datai_x;
Mr_y = datai_x;
Mr_z = datai_x;
FCorrectx = datai_x;
FCorrecty = datai_x;
F_fib = zeros(timeNum(1),y*t*(h-1));
f_Lce = F_fib;
% f_Vce = F_fib;
HillActive = F_fib;
HillPassive = F_fib;
length_per = F_fib;
Vce = F_fib;
Fs_x = zeros(timeNum(1),y*t*h,6);
Fs_y = Fs_x;
Fs_z = Fs_x;
length_s = zeros(timeNum(1),seNum(1));
lengthen_s = length_s;
springF = length_s;
muscleFiberF = length_s;
muscleFiberF_x = datai_x;
muscleFiberF_y = datai_x;
muscleFiberF_z = datai_x;