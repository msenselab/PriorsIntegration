%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Time reproduction to test global priors
% the target intervals are the ramdomly selected elements from the group
% three groups of time intervals the short (400-800ms), long (1200-1600ms),
% and the mixed group(400-800ms, 1200-1600ms). 
% the short and long ranges were separated by 400 ms.
% The Third group is a mixed group of privious two groups. 
% feedback: 1-5, respect to the reproduction error
% left is associated with short and right with long range, or otherwise
% In Exp1, the intervals from short-range were spaced logarithmically between 400 and 800 ms, 
% with a mean of 583 ms, a standard deviation (SD) of 142.2 ms;  the intervals from long-range 
% were spaced logarithmically between 1200 and 2400 ms, with a mean of 1752ms and a standard 
% deviation (SD) of 427.0 ms; distribution of intervals for the mix session had a mean of 1166ms
% and a standard deviation (SD) of 665.5ms. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Exp3_main()

try
    clc
    close all;
    clearvars;
    
    k = 1.2;
    
    [intervalSet1, deltaY] = generateRange(400, 800, 4); %0.3333    0.3964    0.4714    0.5606    0.6667
    intervalSet1 = intervalSet1/(k*1000);

    [intervalSet2] = generateRangeDY(1200, deltaY, 4); %1.4400    1.7125    2.0365    2.4218    2.8800
    intervalSet2 = intervalSet2 * k /1000;
    mean = 1; % classify long and short interval 
    disp(strcat('diff:', num2str(intervalSet2(1)-intervalSet1(end))))  %diff:0.4
    disp( strcat('deltaY:', num2str(deltaY)))  %deltaY:0.17329
    
    KbName('UnifyKeyNames');
       
      kbRP = CInput('m', [1], {'DownArrow'});
    pracBlockNum = 1;
    
    ninterval = (length(intervalSet1))*2; %10
    exp = CExp(4, [3 ninterval], 'blockRepetition', 4); % 480 trials  takes 65 mins %cond group(3);  cond 2 target durations (10) ;
    
    nTrlsBlk = 2*ninterval; %20 trials each block

    intervals = [intervalSet1, intervalSet1; %short interval group
       
    intervalSet2, intervalSet2; %long interval group
                 intervalSet1, intervalSet2];%mixed interval group
            
    % enquire subject information
    subInfo = exp.subjectInfoBox();
    if(strcmp(exp.sName, 'cancelled'))
        disp('This experiment has been cancelled!')
        return;
    end
    %order  = randperm(2); % randromly set the order(shor or long group firstly)
    %set the order
    exp.seq = sortrows(exp.seq, 1);  %sort reward type
    
    if(subInfo.sFirst == 1)
        order = [1, 2, 3];  %short group at first 
    else
        order = [2, 1, 3];  %long group at first 
    end
    for i = 1 : length(order)
        exp.seq(exp.maxTrls/3*(i-1)+1:i*exp.maxTrls/3, 1)  = order(i);
    end
    
    %random generate the target duration
    for i = 1 : exp.maxTrls
          exp.seq(i,3) =  intervals(exp.seq(i,1),exp.seq(i,2));
    end 
    
    
    %a = exp.seq(1 : exp.maxTrls/3, 3); %mean 0.5828 std  0.1422
    %b = exp.seq(exp.maxTrls/3 + 1 : 2*exp.maxTrls/3 , 3);  % mean 1.7485 std  0.4265
    %c = exp.seq(2*exp.maxTrls/3 + 1 : exp.maxTrls , 3);  %mean 1.1656    0.6655
   
    
    %add practice trials
    exp.seq = [exp.seq(1:nTrlsBlk*pracBlockNum,:); exp.seq];
    exp.maxTrls = exp.maxTrls + nTrlsBlk*pracBlockNum;
    
    v = CDisplay('bgColor',48,'fontSize',18,'monitorSize',22, 'fullWindow', 1,'skipSync',1);
    
    infoText = init_text;
    HideCursor;
    %ShowCursor;
    
    % set time reproduction parameters
    % para.nDurations = [0.5  0.8   1.1    1.4    1.7];
    % para.nDurations = durations: [0.5  0.8    1.1    1.4    1.7]; 
    para.vSize = 150; % size is 2 times as large as vM stimuli
    para.fColor = [187,187, 187]; % foreground color
    para.fColorW = [252];
    
    
    %the luminance of the red and green are set to be same.
    para.green = [0, 134, 0];   %need to be adjust for different monitor
    para.red = [192 ,  0, 0];
    para.grey = [128, 128, 128];
    para.circlecolor = {para.red, para.green, para.fColor};
    para.CirleSize = 15;
    para.vFeedDivid = 150;
    para.vFeedd = 15;
    para.xyFeedbackArray = [-2,0; -1, 0; 0, 0; 1, 0; 2,  0] * 1.5;
    para.fbRange = [-100, -0.3; -0.3, -0.05; -0.05, 0.05; 0.05, 0.3; 0.3, 100]; %feedback range with respect to the reproduction error
    vFixation = v.createShape('circle',1, 1,'fill',1, 'color', [250,250,250]);
    vdispPos = [-6, 0; 6, 0; 0, 0]; % atan(12/57)*180/pi()
     
    dotSize = para.vSize/para.vFeedDivid;
    vGreenDisk = v.createShape('circle', dotSize, dotSize, 'color', para.green);
    vRedDisk = v.createShape('circle', dotSize, dotSize, 'color', para.red);    
    vDiskFrame = v.createShape('circle',dotSize, dotSize, 'color', para.fColorW,'fill',0);
    itemSize = 5*dotSize;
    vGreyCircle = v.createShape('circle', itemSize, itemSize,'fill',1, 'color', para.grey);
 
    para.vFullFrames = [vDiskFrame, vDiskFrame, vDiskFrame, vDiskFrame, vDiskFrame];
    para.vFullDisks = [vRedDisk, vGreenDisk, vGreenDisk, vGreenDisk, vRedDisk];
    rectDur = CenterRectOnPoint([0, 0, para.vSize, para.vSize], v.cx, v.cy);

    % display instructions and wait for keypress
    if pracBlockNum >= 1
        pracTex = '\n There are one practice block,';
    else 
        pracTex = ['\n There are ' int2str(pracBlockNum) ' practice blocks,'];
    end 
        
    blockText = [pracTex ' and ' int2str(exp.maxTrls/nTrlsBlk) ' blocks in total.\n\n'];
    v.dispText([infoText.instruction, blockText, infoText.startblock]);

    
    kbRP.wait;
    for trialID = 1:exp.maxTrls
        if pracBlockNum >= 1 && trialID == 1
            v.dispText('Please press any key to start practice block\n');
            kbRP.wait;
        elseif mod((trialID- pracBlockNum*nTrlsBlk-1), nTrlsBlk) == 0
            currentBlockNo = floor((trialID-pracBlockNum*nTrlsBlk)/nTrlsBlk)+1;
            if (currentBlockNo == 9 || currentBlockNo == 17)
                v.dispText(['This part is finished \nPlease take a longer break\n '...
                    'When you are ready press any key to start this block']);
            	kbRP.wait;
            end
            v.dispText(['Block ' num2str(currentBlockNo) '\nPlease press any key to start']);
            kbRP.wait;
        end
        cond = exp.getCondition; %get condition array from exp configuration
        loc = 1;
        if(subInfo.sLoc == 1) 
            loc = (cond(3) > mean)+1; %left is associated with short and right with long range
        else 
            loc = (cond(3) < mean)+1; %left is associated with short and right with long range
        end
        
        phyDuration = 0;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%1. Display a time interval for production-reproduction task
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        v.dispFixation(20);
        WaitSecs(0.5 + rand(1)*0.5); % 500 - 1000 ms 
        v.dispItems(vdispPos(loc,:), vGreyCircle, [itemSize, itemSize]);
        vInitTime = GetSecs;
        WaitSecs(cond(3) - 0.003);
        v.dispText(' ');
        phyDuration = GetSecs - vInitTime;      %visual duration
        v.flip; % blank screen for 1000 ms
        WaitSecs(1);
        
        %%%%%%%%%%%%%%%%%%%%%
        % 2.  reproduction
        %%%%%%%%%%%%%%%%%%%%%
        v.dispText(infoText.reproduction);
        [key, keyInitTime] = kbRP.response; 

        Screen('FillOval', v.wnd,  para.grey , rectDur);
        [vbl, vInitTime] =Screen('Flip', v.wnd);
        keyReleaseTime = kbRP.keyRelease;
        [vbl, vStopTime] =  Screen('Flip', v.wnd);
        repDuration = keyReleaseTime - keyInitTime; % key reproduction
        
        %phyDuration = vStopTime - vInitTime; % visual reproduction
        if kbRP.wantStop
            v.dispText(infoText.stopInfo);
            disp(infoText.stopInfo);
            break;
        end
        % present a feedback display
        feedbackDisplay = para.vFullFrames;
        delta = (repDuration - phyDuration)/phyDuration;
        % find the range of the error
        % para.fbRange = [-100, -0.3; -0.3, -0.05; -0.05, 0.05; 0.05, 0.3; 0.3, 100];
        cIdx = para.fbRange > delta; % column index of left and right boundary
        idx = find(xor(cIdx(:,1),cIdx(:,2)));
        feedbackDisplay(idx(1)) = para.vFullDisks(idx(1));
        
        %WaitSecs(0.25); % wait 250 ms
        v.dispItems(para.xyFeedbackArray, feedbackDisplay,[para.vSize/para.vFeedDivid para.vSize/para.vFeedDivid]); % draw texture
        if abs(delta) > 0.3
            WaitSecs(1.500);
        else
            WaitSecs(0.500); % display the feedback for 500 ms
        end
        
        v.flip;
        WaitSecs(1);
        
        %save responses
        exp.setResp([phyDuration, repDuration]);
        %disp([cond, repDuration]);
         
     end
     if kbRP.wantStop
         v.dispText(infoText.stopInfo);
         disp(infoText.stopInfo);
         Screen('CloseAll');
         return;
     end
    
    %closing the experiment
    exp.saveData;   %save data
    
    v.dispText(infoText.thankyou);
    WaitSecs(1);
    kbRP.wait;
    ShowCursor;
    v.close;
    Screen('CloseAll');
     
catch ME
        ShowCursor;
        disp(ME.message);
        disp(ME.stack);
        for iTrl=1:length(ME.stack)
            disp(ME.stack(iTrl).name);
            disp(ME.stack(iTrl).line);
        end
        v.close;
        Screen('CloseAll');
    end % end try/catch
end % end whole colorworkingmemoryscript



function infoText = init_text
    infoText.instruction = ['During the experiment \n Please follow the instruction the experimenter gave.\n ',...
        'During each trial, you will see a grey circle on the screen, and after a while it disappears. Please feel the time interval of circle. '...
        'The duration of the circle displayed on the screen is the interval you need to replicate.'...
        ' You need to keep pressing right button (Mouse) to reproduce the time duration as long as you think the previous stimulus was presented.'...
        ' After releasing the mouse, you will get a feedback on how good you have reproduced the duration of the stimulus.'...
        ' You can start to reproduce after a downarrow is represented.\n\n'];
	infoText.startblock = ' When you are ready, press any key to start';
    infoText.question = '?';
    infoText.stopInfo = 'Stop in runing experiment...';
    infoText.thankyou = 'The experiment is finished! \n Thank you very much!';
    infoText.redTarget = ['<-'];
    infoText.greenTarget = ['->'];
    infoText.reproduction = ['|\nV'];
end 




