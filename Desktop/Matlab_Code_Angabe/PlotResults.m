function PlotResults()
addpath(genpath('FreeField Results'));

%choose which results to plot, eithe hrtf, freefield or both
result = 'both';

files=dir('/home/bobby/Documents/MSc/MatlabCode/Desktop/Matlab_Code_Angabe/FreeField Results/*.txt');
rod = [1 4 3 2 5];
rod =[rod rod+5 rod+10 rod+15 rod+20 rod+25];

for k=1:length(rod)
    T=load(files(rod(k)).name);
    C = cellstr(strsplit(files(rod(k)).name(1:end-4),'_'));
    mat=vec2mat(T(:,6),10);
    switch result
        case 'freefield'
            h=bar(0:30:180,mat(:,1:7));
            grid on
            l = cell(1,7);
            l{1}='LSFIB'; l{2}='DSB'; l{3}='MVDR Bessel'; l{4}='MVDR Sinc'; l{5}='Kuklasinski Freefield'; l{6}='Kuklasinski';
            l{7}='Unprocessed';
        case 'hrtf'
            h=bar(0:30:180,mat(:,8:10));
            grid on
            l = cell(1,3);
            l{1}='LSFIB\_HRTF'; l{2}='DSB\_HRTF'; l{3}='MVDR\_HRTF';
        case 'both'
            h=bar(0:30:180,mat);
            grid on
            l = cell(1,10);
            l{1}='LSFIB'; l{2}='DSB'; l{3}='MVDR Bessel'; l{4}='MVDR Sinc'; l{5}='Kuklasinski Freefield'; l{6}='Kuklasinski';
            l{7}='Unprocessed'; l{8}='LSFIB\_HRTF'; l{9}='DSB\_HRTF'; l{10}='MVDR\_HRTF';
    end
    
    legend(h,l);
    xlabel('Look Direction');
    ylabel('ASR');
    if (strcmp(C(3),'0') && strcmp(C(4),'0'))
        title(strcat(C(1),'-',C(2),' No Interferer No Noise'));
    elseif (strcmp(C(3),'0') && strcmp(C(4),'1'))
        title(strcat(C(1),'-',C(2),' No Interferer White Noise with 30dB SNR'));
    elseif (strcmp(C(3),'45') && strcmp(C(4),'1'))
        title(strcat(C(1),'-',C(2),' Interferer at ',C(3),' White Noise with 30dB SNR'));
    else
        title(strcat(C(1),'-',C(2),' Interferer at ',C(3),' No Noise'));
    end
    pause
end
end