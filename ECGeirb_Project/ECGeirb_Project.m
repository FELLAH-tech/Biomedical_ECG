%% Main script
% This file computes a simple analysis of an ecg signal. You can use it to test the different processing methods. 
% This first version will plot the temporal signal, compute its cardiac rythma and display the different P, Q, R, S, T points for a specific segment.  
%---------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------Introducion------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------

clear; close all; clc;

%---------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------Data visualization-----------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------

%-------------------------Extract data-------------------------------------
        
        %normal1 = load('C:\Users\hicha\Desktop\ECGeirb_Project\data\ecg_VF.mat');
        normal1 = load('C:\Users\hicha\Desktop\ECGeirb_Project\data\ecg_normal_1.mat');
        x = -normal1.ecg;
        %f_s = load('C:\Users\hicha\Desktop\ECGeirb_Project\data\ecg_VF.mat');
        f_s = load('C:\Users\hicha\Desktop\ECGeirb_Project\data\ecg_normal_1.mat');
        Fs = f_s.Fs;
        
        N_fft = 512;
        window_duration = 4;
        N = window_duration*Fs;
        d = 10;
        w1 = hamming(N);
        
        
%---------------------------Spectrogram------------------------------------
        
        [X1, f1, t1] = stft(x, w1, d, N_fft, Fs);
        S1 = spectro(X1);
        
        imagesc(t1,f1, 10*log10(S1));
        h = colorbar;
        ylabel(h, 'Power/Frequency(dB/Hz)')
        
        xlabel('time(s)');
        ylabel('Frequency(Hz)');
        set(gca,'YDir','normal');
        title('the ecg-normal-1 spectrogram')
        %title('the ecg-VF spectrogram')
        
%--------------------------------------------------------------------------------------------------------------------------------      
%--------------------------------------------PQRST complex Detection-------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
 
    %Data:
        t = 1:1:5;
        Ts = 1/Fs;
        %Low pass
            B = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
            A = [1 -2 1];
        %High pass
            D = [-1 zeros(1,14) 32 -32 zeros(1,14) 1];
            C = [1 -1];
        %Differentiated filter
            E = [1 2 0 -2 -1];
            
    %Filters:
        %Low pass:
            Low_pass = filter(B,A, [1 zeros(1,12)]);
            pass_low = conv(x, Low_pass);
            
        %High pass:
            High_pass = filter(D,C,[1 zeros(1,32)]);
            
        %Band pass:
            pass_band = conv(pass_low, High_pass);
            
        %Differential:
            differential = filter(E, 8*Ts, [1 zeros(1,4)]);
            pass_dif = conv(pass_band, differential);
            
        %squared:
            s_sq = abs(pass_dif).^2 ;
            
        %Normalization
            s_sq_n = s_sq/max(abs(s_sq));
            
        %Moving-window integration:
            N = 30;
            mwi = 1/N*ones(1,N);
            s_mwi = conv(mwi,s_sq);
            s_mwi_n = s_mwi/max(abs(s_mwi));
            
            thresholding = mean(s_mwi_n);
           
            thresh_dec = s_mwi_n > thresholding;
            thresh = decalage(thresh_dec,21,x); % group delay
            
            
            %R wave detection algorithm 
            Max= zeros(1,length(x));
            i = 1;
            while( i <= length(x))
                if(thresh(i)==1) % to locate for each R wave the interval of ones
                    c = 0;
                    while (thresh(i)==1)           
                        if (x(i) > x(i+1) && c == 0 )
                            Max(i) = x(i);
                            c = 1; % to never apply the loop after finding R location
                        elseif (c == 1)
                            Max(i) = 0;
                        end
                        i=i+1;
                    end
                end
                i=i+1;
            end
            
            %Q and S waves detection algorithm
            Q_pos = zeros(1,length(x));
            for i = 1:length(Max)
                if(Max(i)~=0)
                    j = i;
                    if (j-1 > 0)
                        while(x(j-1)<x(j)) % the first minimun before R wave
                            j = j-1;
                        end
                        Q_pos(j) = x(j);
                    end
                end
            end
            
            S_pos = zeros(1,length(x));
            for i = 1:length(Max)
                if(Max(i)~=0)
                    j = i;
                    while(x(j+1)<x(j)) %the first minimun after R wave
                        j = j+1;
                    end
                    S_pos(j) = x(j);
                end
            end
            
            %P and T wave detection algorithm
                %filters
                    G_1 = filter([1 0 0 0 0 0 -1],1,[1 zeros(1,8)]);
                    G1 = conv(G_1,x);
                    G_2 = filter([1 0 0 0 0 0 0 0 -1],[1 -1],[1 zeros(1,10)]);
                    G2 = conv(G_2,G1);
            
                %R-R interval
                    R_R = zeros(1,length(Max)); %to adapt it with the vector "Max"
                    R_R2 = []; % used for bpm
                    i = 1;
                    while( i < length(Max))
                        if( Max(i)~=0) % R location
                            c = 1; % to compute the R-R interval
                            p = i; % the R location index used for the vector "R_R
                            while (Max(i+1)==0)
                                c = c + 1;
                                i = i+1;
                                if (i == length(Max))%to Avoid exceed the length 
                                    break;
                                end
                            end
                            R_R(p) = c;
                            R_R2 = [R_R2 c];
                        end
                        i = i+1;
                    end
                    
                %R-T interval
                    R_T = 0.7*R_R;
                
                %P and T detection 
                    Max_dec = decalage(Max,-6,G2); % group delay
                    R_T_dec = decalage(R_T,-6,G2);
                    R_R_dec = decalage(R_R,-6,G2);
                    T_pos_dec = zeros(1,length(x));
                    P_pos_dec = zeros(1,length(x));
                    i = 1;
                    while( i < length(Max_dec))
                        if( Max_dec(i)~=0)
                                %T detection
                                    verif = 0;
                                    for k = i+1:floor(R_T_dec(i)+i+1)
                                       if ( k < length(S_pos)) %to Avoid exceed the length 
                                           if(S_pos(k)~=0 || verif == 1 ) % to avoid having a T before S
                                                verif = 1;
                                                if (G2(k-1)>0 && G2(k+1)<0) % interval where G2 passed by zero
                                                    T_pos_dec(k)= x(k);
                                                end
                                           end
                                       end
                                    end
                                    max_T = T_pos_dec(i+1);
                                    
                                %searching for the hightest value in ECG 
                                    for k = i+1:floor(R_T_dec(i)+i+1)
                                        if ( k < length(T_pos_dec))
                                            if (max_T < T_pos_dec(k+1))
                                                T_pos_dec(k)= 0;
                                                max_T = T_pos_dec(k+1);
                                            elseif (max_T > T_pos_dec(k+1))
                                                T_pos_dec(k+1)= 0;
                                            end
                                        end
                                    end
                                    if ( floor(R_T_dec(i)+i+2) < length(P_pos_dec))
                                        max_P = P_pos_dec(floor(R_T_dec(i)+i+2));
                                    end
                                %P detection
                                    pos_max = 1; %to avoid an error if the pos_max is not found
                                    for k = floor(R_T_dec(i)+i+2):R_R_dec(i)+i % the 0.3 R-R interval remaining
                                       if ( k < length(P_pos_dec))
                                           if (max_P < x(k+1))
                                              max_P = x(k+1);
                                              pos_max = k+1;
                                           end
                                           if (Q_pos(k)~=0)
                                                break;
                                           end
                                       end
                                    end
                                    if (pos_max~=0) %to avoid an error if the P wave is not detected
                                    P_pos_dec(pos_max)=max_P;
                                    end
                        end
                        i=i+1;
                    end
                %The first P wave detection 
                for k = 1:i
                    if ( k < length(P_pos_dec))
                       if (max_P < x(k+1))
                           max_P = x(k+1);
                           pos_max = k+1;
                       end
                       if (Q_pos(k)~=0)
                           break;
                       end
                    end
                end
                if (pos_max~=0) %to avoid an error if the P wave is not detected
                P_pos_dec(pos_max)=max_P;
                end
                
                T_pos = decalage(T_pos_dec,6,x);%goup delay
                
                
                %Tachycardia/Bradycardia
                    mean_R_R = sum(R_R2)/(length(R_R2)+1);
                    bps = mean_R_R/Fs;
                    bpm = 60/bps;
                
                %Ectopic beat detection algorithm
                    Ectopic_detected = 0;
                    epsilon = 3*mean_R_R;
                    for i = 1:length(R_R2)-1
                        if(abs(R_R2(i)-R_R2(i+1)) > epsilon)
                            Ectopic_detected = 1;
                        end
                    end
                
                %Atrial fibrillation 
                    ACF = zeros(1,200);
                    for k=1:length(ACF)
                        ACF(k) = autocovariance(k,R_R2,mean_R_R);
                    end
                    
 %---------------------------------------------------------------------------------------------------------------------------
 %-------------------------------------------------Figures-------------------------------------------------------------------
 %---------------------------------------------------------------------------------------------------------------------------
 
        %illustrate locations in the figure
            figure, plot(x), xlim([0 1000]), title('ECG segment characteristic'),grid MINOR;
            figure, plot(x), xlim([0 400]), title('ECG segment characteristic'),grid MINOR;
            xlabel('temps (s)');
            ylabel('Magnitude');
            hold on 
            for i=1:length(x)
                if(S_pos(i)~=0)
                    plot(i,x(i),'*','Color','red');text(i,x(i),' S ','Color','red','FontSize',14);
                end
                
                if(Max(i)~=0)
                    plot(i,x(i),'*','Color','red');text(i,x(i),' R ','Color','red','FontSize',14);
                end
                
                if(Q_pos(i)~=0)
                    plot(i,x(i),'*','Color','red');text(i,x(i),' Q ','Color','red','FontSize',14);
                end
                
                if(T_pos(i)~=0)
                    plot(i,x(i),'*','Color','red');text(i,x(i),' T ','Color','red','FontSize',14);
                end
                
                if(P_pos_dec(i)~=0)
                    plot(i,x(i),'*','Color','red');text(i,x(i),' P ','Color','red','FontSize',14);
                end
            end
            hold off;               
            
        %Bpm comparison
            figure,plot(bpm*ones(1,100));hold on;plot(60*ones(1,100),'g');plot(100*ones(1,100));hold off
            title('The Tachycardia/Bradycardia detection (normal ECG)');
            %title('The Tachycardia/Bradycardia detection (PVC ECG)');
            legend('bpm','bradycardia threshold','Tachycardia threshold');
            ylabel('beats per minute');
            ylim([0 150]);
            xlim([1 50]);
            
        %Ectopic detection
            figure,plot(Ectopic_detected*ones(1,100));
            title('The Ectopic detection (normal ECG)');
            ylabel('Ectopic detected if 1');
            ylim([-1 2]);
            xlim([1 50]);
            
        %Auto-covariance function (ACF)
            figure,plot(ACF);grid on
            ylabel('Auto-covariance Function');
            xlabel('k sequence');
            title('Auto-covariance Function for AF ECG');
            
    
