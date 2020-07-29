function [RR_out] = remove_ect(RR, thresh1)

% [RR_out] = remove_ect(RR, thresh1, thresh2, samp_freq);
%
% Takes the two col RR data (t, x) and removes any RR interval
% <> 80% of the previous (recorded) interval. mean, std and kutosis are
% recorded and a histogram and time series are plotted for each 5
% minute section (overlapping by 4 mins)


if nargin < 2
    thresh1=0.8;
end



RR_t = RR(:,1);
RR_x = RR(:,2);
%%%%%
count =2;
RR_t_new(1)=RR_t(1);
RR_x_new(1)=RR_x(1);
post_ect_flag = 0;
last_real_RR=RR_x(1);
post_pef=0;

for y=2:length(RR_t) % 1: the length of the 5 minute buffer
    %a(y)=abs(RR_x(y)-last_real_RR);
    %b(y)=(1-thresh1)*RR_x(y-1)  ;
    if(post_pef==1) % if post post ectopic beat
        post_ect_flag=0; % set flags to use
        post_pef=0;      % next datum
    end
    
    if post_ect_flag==1   % if the last beat was ectopic skip this loop and
        % use it only for calculating the next RR interval
        post_pef=1;
    end
    
    if post_ect_flag==0 % ignoring post ectopic beats
        % IF the new one is more than SDout SDs outside mean of entire
        %	data reject it --- not very adaptive !!!!
        
        % if(abs(RR_t_raw_buff(y)-mean_of_data)<SDout*std_of_data)
        % if new beat is <>80% of previous beat, ignore it.
        if abs(RR_x(y)-last_real_RR)<=((1-thresh1)*RR_x(y-1))
            RR_x_new(count)=RR_x(y);
            RR_t_new(count)=RR_t(y);
            last_real_RR=RR_x(y);
            count=count+1;
        else
            disp('ectopic')
        end
        if abs(RR_x(y)-last_real_RR)>((1-thresh1)*RR_x(y-1))
            post_ect_flag=1;
        end
        
    end
    
end

RR_out = [RR_t_new' RR_x_new'];

%figure;

%plot(RR_t,RR_x,'+--b');
%hold on
%plot(RR_t_new,RR_x_new,'*r');
%xlabel('sample no'), ylabel('RR');
%%axis([0 inf 30 200])
%%title(['Time series- 5mins of data without ectopics or followers, Chunk: ' num2str(j)]);
%figure
%histfit(RR_x_new,20)
%xlabel('HR'), ylabel('frequency');
%%axis([30 200 0 100])
%%title(['PSD of 5mins of data without ectopics or followers, Chunk: ' num2str(j)]);

%figure;
%hold on
%plot(a,'.--');
%plot(b,'.--r');







