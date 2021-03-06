
%name=1;
ratio=[];
emgpow=[];
name
eegband=['delta';'theta';'alpha';'beta_';'gamma'];
for datas = 1:3
        datas 
        if datas==1
            data=load(strcat('../',int2str(name),'/rest1_RawData.dat'));
            duration=240;
        elseif datas==3
            data=load(strcat('../',int2str(name),'/rest2_RawData.dat'));
            duration=240;
        else
            data=load(strcat('../',int2str(name),'/sport_RawData.dat'));
            duration=540;
        end

       
        c1=data(1,:)-(data(1,:)+data(2,:)+data(3,:)+data(4,:))/4;
        quad=data(8,:);
        heart=-data(9,:);

        ndata(1,:)=c1;
        ndata(2,:)=quad;

        fmin    =   [ 1 , 4 , 8 , 12, 25];
        fmax    =   [ 4 , 8 , 12, 25, 60];
        fres    =   [ 1 ,  1 ,  1 , 1, 1];
        start   =   1;
        over    =   size(ndata,2);
        samplerate = 1000;

        m=[];

        [TD LE] =   cutting(ndata,start,over,samplerate);

        segLen=5*samplerate;

        %%% ECG find RRinterval %%%

            % remove drift %
            heart2=timeseries(heart',(1:length(heart))/1000);
            heart3 = idealfilter(heart2, [1 100], 'pass');
            ECG=heart3.Data';
            %%%

        minPeakH=800;
        minDis=300;
        [~,locs_R] = findpeaks(ECG,'MinPeakHeight',minPeakH,'MINPEAKDISTANCE',minDis);
        RRfreq=[];
        locs_R=[locs_R(3:length(locs_R)-1)];
        RRfreq=[RRfreq 60*samplerate/(locs_R(2)-locs_R(1))];
        for i = 2 : size(locs_R,2)-1
            RRfreq=[RRfreq 60*samplerate*2/(locs_R(i+1)-locs_R(i-1))];
        end
        RRfreq=[RRfreq 60*samplerate/(locs_R(size(locs_R,2))-locs_R(size(locs_R,2)-1))];

            % RR interval average %
            x=10;
            RRave=[];
            for time = 1 : 5*samplerate : duration*samplerate
                RRtmp=[];
                while( x<length(locs_R) && locs_R(x)<=time+5*samplerate )
                    RRtmp=[RRtmp RRfreq(x)];
                    x=x+1;
                end
                RRave=[RRave mean(RRtmp)];
            end
            %%%

        %%%%%%
        
        %%%EEG ratio%%%
        channel=1;
        for band = 1 :length(fmin) 
            band
            pow=[];
            filtwave=[];
            for time = 1 : 5*samplerate : duration*samplerate
                filtwave=mean(tfa_morlet(TD(channel, time : time+segLen),samplerate,fmin(band),fmax(band),fres(band)));
                pow = [pow  mean(filtwave)];
            end
            POW(band,:)=pow;
            save(strcat('../',int2str(name),'/',eegband(band,:)),'POW');
        end
        
        for time = 1 : size(POW,2)
            ratio = [ratio (POW(1,time)+POW(4,time)+POW(5,time))/POW(3,time)];
        end
        
        %%%EMG%%%
        channel=2;
        filtwave=[];
        for time = 1 : 5*samplerate : duration*samplerate
            filtwave=mean(tfa_morlet(TD(channel, time : time+segLen),samplerate, 50 , 150 , 1 ));
            emgpow = [emgpow  mean(filtwave)];
        end
        
        if datas==1
            ave1=mean(ratio);
            ave2=mean(emgpow);
        end
        
        clearvars -except ratio emgpow name ave1 ave2 eegband
end

ratio=ratio/ave1;
emgpow=emgpow/ave2;

save(strcat('../',int2str(name),'/ratio'),'ratio');
save(strcat('../',int2str(name),'/emg150'),'emgpow');

clear
