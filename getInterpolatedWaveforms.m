function [lactate_raw,varargout]=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,show,average_window,realTimeFlag)
%
% waveforms= getInterpolatedWaveforms(varLabels,category)
%
% Generates interpolated waveforms from unevenly sampled data
% in category.
%
% varLabeLs Nx1 cell array of the physiological labels to category of the waveforms to be
% interpolated
%
% category Mx3 cell array
%
% Ts 1x1 sampling interval (in hours) from which to interpolate data
%
% outVarName N+1x1 cell array for the names of the new interpolated data
% First column is the raw lactate values
%
% Output is multiple variables, one for each element of varLabels see
% example below:
%
% %Example 1:
% outVarName={'lact','map','hr','urine','weight','pressor'};
%[lact,map,hr,urine]=getInterpolatedWaveforms(varLabels,category,Ts,outVarName);

NvarName=length(varLabels);
lactate_raw=[];
if(isempty(realTimeFlag))
    %Default =0 will use FILFILT for all filtering operations to have 0
    %phase
    realTimeFlag=0;
end
%Average window in units of hour
if(isempty(average_window))
    average_window=0.5;
end
if(average_window==0)
    b=1;
    nb=1;
else
    nb=round(average_window/Ts); %Filter waveforms with half an hour moving average
    b=ones(nb,1)./nb;
end
xraw=[];
for n=1:NvarName
    
    ind=strcmp(category,varLabels{n});
    if(ind==0)
        %No data for this time series, continue
        eval([outVarName{n} '=[NaN NaN];'])
        continue
    end
    
    x=[tm(ind) val(ind)];
    del=find(isnan(x(:,1))==1);
    if~(isempty(del))
        x(del,:)=[];
    end
    isUrine=1; %Urine requires special attention.
    if(~strcmp(varLabels{n},'URINE'))
        %Remove cases where value is zero for all time series except urine
        %and the ones in the if condition below
        if(~strcmp(varLabels{n},'PRESSOR_TIME_MINUTES') && ~strcmp(varLabels{n},'MECH_VENT_FLAG'))
            del=find(x(:,2)==0);
            if(~isempty(del) && ~strcmp(varLabels{n},'PaCO2'))
                warning(['Time series has zeros values : ' varLabels{n}])
            end
        end
        isUrine=0; %set flag for false, this is not urine series
    end
    [Nx,~]=size(x);
    if(show)
        xraw=x;
    end
    if(strcmp(varLabels{n},'LACTATE'))
        lactate_raw=x; 
    end
    
    if(Nx>1)
        x=sortrows(x,1);
        %Average any points measue at the same time
        diff_tm=diff(x(:,1));
        ind_rep=find(diff_tm==0);
        if(~isempty(ind_rep))
            if(isUrine)
                %In urine case, sum the points together instead of
                %averaging
                for reps=1:length(ind_rep)
                    x(ind_rep+1,2)=x(ind_rep+1,2)+x(ind_rep,2);
                    x(ind_rep,2)=NaN;
                end
                %Remove the repeats
                del=find(isnan(x(:,2))==1);
                x(del,:)=[];
            else
                %For other series, take the approximate average
                %(this will be biased towards the last sample)
                for reps=1:length(ind_rep)
                    x(ind_rep+1,2)=mean(x(ind_rep+1,2)+x(ind_rep,2));
                    x(ind_rep,2)=NaN;
                end
                %Remove the repeats
                del=find(isnan(x(:,2))==1);
                x(del,:)=[];
            end
            
        end
        
        if(isUrine)
            %Do not apply median to urine, calculative cumulative hourly
            x=hourly_cumulative(x);
        else
            x=hourly_median(x);
        end
        
        %Interpolate waveforms
        [Nx,~]=size(x);
        if(Nx==0)
            y=[NaN NaN];
        elseif(Nx==1)
            y=x;
        else
            y=[x(1,1):Ts:x(end,1)]';
            y(:,2)=interp1(x(:,1),x(:,2),y,'linear');
            [Ny,~]=size(y);
            if((Ny+1)> (nb*3))
                %Filter the waveforms through a moving average
                %filtfilt only works for cases where Ny is 3x filter order
                if(nb ~=1)
                    %TODO: Consider using FILTER instead
                    if(realTimeFlag)
                        y(:,2)=filter(b,1,y(:,2));
                    else
                        y(:,2)=filtfilt(b,1,y(:,2));
                    end
                end
            end
        end
    else
        if(Nx==1)
            %Case where there is only one value of x
            y=x;
        else
            %Nx=0
            y=[NaN NaN];
        end
    end
    
    %Set y to the time series being analyzed
    eval([outVarName{n} '=y;'])
    if(show)
        figure
        plot(xraw(:,1),xraw(:,2),'o','LineWidth',2)
        hold on;grid on
        plot(y(:,1),y(:,2),'r','LineWidth',2)
        title(outVarName{n})
        legend('Raw Data','Interpolated Series')
        xlabel('Time (hours')
        ylabel('Measurement Values')
    end
end
for n=1:nargout-1
    eval(['varargout{n}=' outVarName{n} ';'])
end
