function[] = select_5steps(datafolder,subjects,trials)
%Use select5steps to select 5 step range for each subject and trial
addpath(genpath(datafolder))
cd(datafolder)

% start with what we have already and add to that
load('5steps_indices.mat','hsl','start')

for subj = subjects
    disp(subj)
    
    for trial = 1:length(trials)
        disp(trials(trial))
        netshift = 0;
        shift = 10;

        if exist(['p',num2str(subj),'_RawData.mat'],'file')
            load(['p',num2str(subj),'_RawData.mat'],'data')
        else, continue
        end
    
        if isempty(data(trials(trial)).grf)
            continue
        end
        
        f1= data(trials(trial)).grf.force1;
        f2= data(trials(trial)).grf.force2;
        
        %ask if hsl exists or not, if it exists use previous range, if not
        %pop up figure to select heekstrike range
        if isnan(start(subj, trials(trial)))
            figure; 
            plot(f1(1)); hold on; plot(f2{1})
            xlabel ('Force')
            ylabel ('Frames')
            %if isnan(start(subj, trial))
            [idx,~] = ginput(2);
            close
            %idx=both clicks
            idx = round(idx);
            
        else idx(1)=start(subj, trials(trial));
             idx(2)=length(f1(:,1));
        end
        
        % extract ground reaction force
        grfl = f1(idx(1):idx(2),:);
        grfr = f2(idx(1):idx(2),:);

        % Threshold at 20 N
        grflt = grfl;
        grfrt = grfr;
        grflt(grflt(:,3)<20,:) = 0;
        grfrt(grfrt(:,3)<20,:) = 0;  

        idx_mo = round(idx/10);
%pjl and pjr are power
        pjl = [sum(data(trials(trial)).Link_Model_Based.l_ank_power(idx_mo(1):idx_mo(2),:),2)...
               sum(data(trials(trial)).Link_Model_Based.l_kne_power(idx_mo(1):idx_mo(2),:),2)...
               sum(data(trials(trial)).Link_Model_Based.l_hip_power(idx_mo(1):idx_mo(2),:),2)];
        pjr = [sum(data(trials(trial)).Link_Model_Based.r_ank_power(idx_mo(1):idx_mo(2),:),2)...
               sum(data(trials(trial)).Link_Model_Based.r_kne_power(idx_mo(1):idx_mo(2),:),2)...
               sum(data(trials(trial)).Link_Model_Based.r_hip_power(idx_mo(1):idx_mo(2),:),2)];
           
        % extract heelstrikes and toe-offs from grf
        hsl_new = unique(invDynGrid_getHS_TO(grflt, grfrt, 20));
        hsl_mo = unique(round(hsl_new/10) + 1); % heelstrikes for mocop (fs_mocop = .1*fs_grf)
        
        %n= #of heelstrikes, if the pre selected heelstrike exists use that
        %as the start of the 5 step range, if not use the midpoint as the
        %beginning of 5 steps
        n = length(hsl_new);
        if isnan(hsl(subj, trials(trial)))
            netshift=0
        else
            [~,i] = min(abs(hsl(subj, trials(trial))-hsl_new));
                netshift = (i)-(ceil(n/2));
            end

    while shift ~= 0
        close all
       
%         keyboard
        %netshift=0 at the beginning so 5 step start at middle of hsl_new,
        %unless the previously selected range is used, in which case
        %netshift~= 0
        if (netshift) > abs(ceil(n/2)),
           hsl_mid = hsl_new(1:5);
        elseif (netshift) < (ceil(n/2)*(-1)),
            hsl_mid=hsl_new(n-5:n);
        else, hsl_mid = hsl_new((ceil(n/2)+netshift):(ceil(n/2)+5+netshift)); 
        end
        
        hsl_mo_mid = round(hsl_mid/10) + 1;
    
        h = figure;
        
        %% left
        %plot GRF of entire stance phase and GRF of 5 steps
        subplot(241);
        plot(pjl); hold on
        
        plot([hsl_mo_mid(1) hsl_mo_mid(1)], [min(min(pjl)) max(max(pjl))],'k','linewidth',2)
        plot([hsl_mo_mid(end) hsl_mo_mid(end)], [min(min(pjl)) max(max(pjl))],'k','linewidth',2)
        ylabel ('Power(W)')
        
        subplot(242);
        plotperstride(pjl, hsl_mo); hold on
        plotperstride(pjl, hsl_mo_mid)
        ylabel ('Power(W)')
        
        subplot(243);
        plotperstride(pjl, hsl_mo_mid)
        ylabel ('Power(W)')
        
        subplot(244);
        plotperstride(grfl, hsl_mid)
        ylabel ('Force(N)')
        
        %% right
        subplot(245);
        plot(pjr); hold on
        
        plot([hsl_mo_mid(1) hsl_mo_mid(1)], [min(min(pjr)) max(max(pjr))],'k','linewidth',2)
        plot([hsl_mo_mid(end) hsl_mo_mid(end)], [min(min(pjr)) max(max(pjr))],'k','linewidth',2)
        ylabel ('Power(W)')
        
        subplot(246);
        plotperstride(pjr, hsl_mo); hold on
        plotperstride(pjr, hsl_mo_mid)
        ylabel ('Power(W)')
        
        subplot(247);
        plotperstride(pjr, hsl_mo_mid)
        ylabel ('Power(W)')
        
        subplot(248);
        plotperstride(grfr, hsl_mid)
        ylabel ('Force(N)')
        
        for i=1:8
            subplot (2,4,i)
            xlabel ('frames')
            axis tight
            
        end
        
        set(h,'Units','normalized', 'position', [0 .5 1 .5])
        
        %Ask user if they want to shift, the number entered shifts them
        %that many heelstrikes forward or backward from the current 5 step
        %range (not the first 5 step range)
        shift = input('Want to shift? [0 = no shift, + = forward shift, - = backward shift]');
        close(h)
        
        netshift = netshift + shift;

    end

    hsl(subj,trial) = hsl_mid(1);
    start(subj,trial) = idx(1);
    end
end
    
save('5steps_indices.mat','hsl','start')

end