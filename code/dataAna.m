%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name: dataAna.m
%%
%% Author: Fiona Zhu
%% Date: Jan 2019
%% Data analyse code for the Time reproduction task mixed with WM task
%% to test serial dependence between trials and within a trial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note that function 'sgtitle' is supported in R2018b. please comment the
%code using 'sgtitle', if it is not runable on you computer.

function dataAna

    clear all; close all;  clc; 
    
    allArray_S1= readtable('../data/AllData_study1.csv'); 
    allArray_S1 = allArray_S1(allArray_S1.valid == 1, :);
   
    %allArray_S2= readtable('../data/AllData_study2.csv');
    allArray_S3= readtable('../data/AllData_study3.csv');
    allArray_S3 = allArray_S3(allArray_S3.valid == 1, :);
    maxX = max(allArray_S1.firstDur);
    xlist = -maxX: 0.05 :maxX;
    %xlist = -maxX-4: 0.05 :maxX+4; 
    
  

    
    par_seq1_raw =  SequentialEffectfitRawDat(allArray_S1, 'S1');
    %par_seq2_raw =  SequentialEffectfitRawDat(allArray_S2, 'S2');
    par_seq3_raw =  SequentialEffectfitRawDat(allArray_S3, 'S3');
    
    par_con1_raw = ShowContrastEffectRawDat(allArray_S1, 'S1');
    %par_con2_raw = ShowContrastEffectRawDat(allArray_S2, 'S2');
    par_con3_raw = ShowContrastEffectRawDat(allArray_S3, 'S3');

    plotfitDOG(xlist, par_seq1_raw, par_con1_raw, par_seq3_raw, par_con3_raw,'(Raw Data)');
    
    
%      % uncommment the following codes to run fitting on rund data
%      %processing data 
%      allArray_S1 = datprocessing(allArray_S1);   %round 
%      %allArray_S2 = datprocessing(allArray_S2);   %round
%      
%      plotRep_repError(allArray_S1, 'S1'); 
%     %plotRep_repError(allArray_S2, 'S2'); 
% 
%     %showAssimilationEffect(allArray_S1, 'S1'); 
%     %showAssimilationEffect(allArray_S2, 'S2'); 
%     
%     par_seq1 = plotSequentialEffect(allArray_S1, 'S1'); 
%     %par_seq2 = plotSequentialEffect(allArray_S2, 'S2');
%     
%     par_con1 = ShowContrastEffect(allArray_S1, 'S1');
%     %par_con2 = ShowContrastEffect(allArray_S2, 'S2');
% 
%   
%     plotall(xlist, par_seq1, par_con1, 'S1');
%     %plotall(xlist, par_seq2, par_con2, 'S2');
%     
%     %plotfitDOG(xlist, par_seq1, par_con1,par_seq2, par_con2, '(Round)');
%    
    
    
end


function plotfitDOG(xlist, par_seq1, par_con1,par_seq2, par_con2, filename)
 figure; 
    hold on;
    plot(xlist, gaussmf(xlist,[par_seq1(2) par_seq1(1)])-gaussmf(xlist,...
                [par_seq1(2) -par_seq1(1)]), 'b', 'LineWidth', 3);
    plot(xlist, gaussmf(xlist,[par_con1(2) par_con1(1)])-gaussmf(xlist,...
                [par_con1(2) -par_con1(1)]), 'r', 'LineWidth', 3);
   
    plot(xlist, gaussmf(xlist,[par_seq2(2) par_seq2(1)])-gaussmf(xlist,...
                [par_seq2(2) -par_seq2(1)]), 'b--', 'LineWidth', 3);
    plot(xlist, gaussmf(xlist,[par_con2(2) par_con2(1)])-gaussmf(xlist,...
                [par_con2(2) -par_con2(1)]), 'r--', 'LineWidth', 3);
   
    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('Relative deviation to previous stimulus');
    ylabel('Repr. Error (secs)');
    legend('Sequential effect(cross trials) (S1)', 'Contrast effect(within trial) (S1)',...
        'Sequential effect(cross trials) (S2)', 'Contrast effect(within trial) (S2)','Location','northeast');
    title(strcat('Sequentual effect and contrast effect (S1+S2)',{' '},filename));
    saveas(gcf,strcat('../figures/ALL', filename, '.png'));

end


function plotall(xlist, par_seq, par_con, filename)
   %% visualization of sequential effect and contrast effect 
    figure; 
    hold on;

    plot(xlist, gaussmf(xlist,[par_seq(2) par_seq(1)])-gaussmf(xlist,...
                [par_seq(2) -par_seq(1)]), 'b', 'LineWidth', 3);
    
    plot(xlist, gaussmf(xlist,[par_con(2) par_con(1)])-gaussmf(xlist,...
                [par_con(2) -par_con(1)]), 'r', 'LineWidth', 3);
   
    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('Duration (secs)');
    ylabel('Repr. Error (secs)');
    legend('sequential effect(cross trials)', 'contrast effect(within trial)','Location','northeast');
    title(strcat('Sequentual effect and contrast effect',{' '}, filename));
    saveas(gcf,strcat('../figures/Sequentual+Contrast', filename,'.png'));
end 


function allArray = datprocessing(allArray)
    %round physical presented duration (100 ms)
    allArray.dDurRound = round(allArray.dDur *10)/10;  %round dDur 
    allArray.nttDur_Round = round(allArray.nttDur * 10)/10; % round nttDur 
    allArray.secondDurRound = round(allArray.secondDur *10)/10; 
    allArray.firstDurRound = round(allArray.firstDur *10)/10; 
    allArray.duration = round(allArray.targetDur *10)/10; 

end

function plotRep_repError(allArray, filename)
    allArray.repError = allArray.RP - allArray.targetDur; 
    avg = grpstats(allArray, {'NSub','duration','target'},{'mean','sem'}, ...
        'DataVars', {'repError','RP'});
    nsub = length(unique(avg.NSub));
    lineColor = ['r', 'g', 'b'];
    
    %% visualization mean_repError
    % individual participant   
    figure;
    for i = 1:nsub 
        subplot(4,4,i);
        hold on;
        idx = avg.NSub == i & avg.target == 1; %red
        plot(avg.duration(idx), avg.mean_repError(idx),'r-o');
        idx = avg.NSub == i & avg.target == 2; %green
        plot(avg.duration(idx), avg.mean_repError(idx),'g-d');
        xlabel('Target Durations (secs)');
        ylabel('Reproduction Error (secs)');
        legend('Red', 'Green','Location','northeast');
    end
    hold off;
    sgtitle(strcat('Error of each participant',{' '}, filename));
    saveas(gcf, strcat('../figures/Error_NSub', filename, '.png'));



    %% visualization mean_repError
    % grand average
    all_avg = grpstats(avg, {'target','duration'},{'mean','sem'}, ...
        'DataVars', {'mean_repError','mean_RP'});
    figure; 
    hold on;
    plot(all_avg.duration(all_avg.target == 1), all_avg.mean_mean_repError(all_avg.target == 1),'r-o');  
    plot(all_avg.duration(all_avg.target == 2), all_avg.mean_mean_repError(all_avg.target == 2),'g-d');
    hold off;
    xlabel('Target Durations (secs)');
    ylabel('Reproduction Error (secs)');
    legend('Red', 'Green','Location','northeast');
    set(gca,'fontsize',8)
    title(strcat('Average Error',{' '}, filename));
    saveas(gcf,strcat('../figures/Error_Average', filename,'.png'));
    
    
  %% visualization of reproducation 
  % individual participant   
    figure;
   
    for i = 1:nsub 
        fig = [];
        subplot(4,4,i);
        hold on;
        for  k = 1:2
            idx = avg.NSub == i & avg.target == k; %red    
            x1 = avg.duration(idx);
            y1 = avg.mean_RP(idx);

            pfit = polyfit(x1, y1, 1);  
            x2 = linspace(min(x1)-1,max(x1)+1);  
            y2 = polyval(pfit,x2);  
            if (k == 1)
                p1 = plot(x1, y1, 'r.'); 
                p2 = plot(x2, y2, 'r', 'LineWidth', 2);
            else 
                p1 = plot(x1, y1, 'gd'); 
                 p2 = plot(x2, y2, 'g', 'LineWidth', 2);
            end          
            fig = [fig  p1 p2];
        end 
        xlabel('Target Durations (secs)');
        ylabel('Reproduction (secs)');
        legend(fig, 'Red Target', 'Red(fit)', 'Green Target', 'Green(fit)','Location','northeast');
    end
    hold off;
    set(gca,'fontsize',8)
    sgtitle(strcat('Reproduction',{' '}, filename));
    saveas(gcf,strcat('../figures/reproducation_NSub', filename,'.png'));
    
    %% visualization of reproducation 
    % grand average
    figure; 
    fig = [];
    hold on;
    for  k = 1:2
        x1 = all_avg.duration(all_avg.target == k);
        y1 = all_avg.mean_mean_RP(all_avg.target == k);
        if (k == 1)
             p1 = plot(x1, y1, 'r.'); 
        else 
             p1 = plot(x1, y1, 'gd'); 
        end
        pfit= polyfit(x1, y1, 1);  
        x2 = linspace(min(x1)-1,max(x1)+1);  
        y2= polyval(pfit,x2); 
        p2 = plot(x2, y2, lineColor(k), 'LineWidth', 2);
        fig = [fig  p1 p2];
    end
    plot([0,2.4],[0,2.4],'--'); % diag line
    hold off;
    xlabel('Target Durations (secs)');
    ylabel('Reproduction (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'Location','northeast');
    title(strcat('Average Reproduction',{' '}, filename));
    set(gca,'fontsize',8)
    saveas(gcf,strcat('../figures/reproducation', filename,'.png'));
    
end


%% visualization to show assimilation effect (D2)
function showAssimilationEffect(allArray, filename)
    
    avg_D2 = grpstats(allArray, {'NSub','secondDurRound','target'},{'mean','sem'}, ...
        'DataVars', {'repError','RP'});
    nsub_D2 = length(unique(avg_D2.NSub));
    lineColor = ['r', 'g', 'b'];

    figure;
    for j = 1:nsub_D2 
        fig = [];
        subplot(4,4,j);
        hold on;
        for  k = 1:2
            idx = avg_D2.NSub == j & avg_D2.target == k; 
            x1 = avg_D2.secondDurRound(idx);
            y1 = avg_D2.mean_RP(idx);
            p=polyfit(x1, y1, 1);  
            x2 =linspace(min(x1)-1,max(x1)+1);  
            y2=polyval(p,x2);  
            if (k == 1)
                 p1 = plot(x1, y1, 'r.'); 
            else 
                 p1 = plot(x1, y1, 'gd'); 
            end
            p2 = plot(x2, y2, lineColor(k), 'LineWidth', 2);
            fig = [fig  p1 p2];
        end 
        plot([0,2.4],[0,2.4],'--'); % diag line
        xlabel('Second Duration (secs)');
        ylabel('Reproduction (secs)');
        legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'Location','southeast');
    end
    hold off;
    sgtitle(strcat('Assimilation effect', {' '}, filename ));
    saveas(gcf,strcat('../figures/assimilation_NSub', filename,'.png'));

    
     %% grand average
    all_avg_D2 = grpstats(avg_D2, {'target','secondDurRound'},{'mean','sem'}, ...
        'DataVars', {'mean_repError','mean_RP'});
    figure; 
    fig = [];
    hold on;
    for  k = 1:2
        x1 = all_avg_D2.secondDurRound(all_avg_D2.target == k);
        y1 = all_avg_D2.mean_mean_RP(all_avg_D2.target == k);
        p = polyfit(x1, y1, 1);  
        x2 = linspace(min(x1)-1,max(x1)+1);  
        y2 = polyval(p,x2);  
        if (k == 1)
             p1 = plot(x1, y1, 'r.'); 
        else 
             p1 = plot(x1, y1, 'gd'); 
        end
       
        p2 = plot(x2,y2, lineColor(k), 'LineWidth', 2);   
        fig = [fig  p1 p2];
    end
    plot([0,2.4],[0,2.4],'--'); % diag line
    hold off;
    xlabel('Sencond duration (secs)');
    ylabel('Reproduction (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'Location','southeast');
    title(strcat('Assimilation effect',{' '}, filename));
    saveas(gcf, strcat('../figures/assimilation', filename,'.png'));
    
end 


%% visualization to show sequential effect for raw data
function par_seq = SequentialEffectfitRawDat(allArray, filename)
    lineColor = ['r', 'g', 'b'];

    nsub = length(unique(allArray.NSub));
    xlist = min(allArray.dDur): 0.05 :max(allArray.dDur);
    figure;

    for j = 1:nsub 
        fig = [];
        subplot(4,4,j);
        hold on;
        for k = 1:2     
            idx = allArray.NSub == j & allArray.target == k; 
            x = allArray.dDur(idx); 
            y = allArray.repError(idx); 
           
            if (k == 1)
                 p1 = plot(x, y, 'r.'); 
            else 
                 p1 = plot(x, y, 'g.'); 
            end
            
            par1 = fitDOG(x, y);
            p2 = plot(xlist, gaussmf(xlist,[par1(2) par1(1)])-gaussmf(xlist,...
                [par1(2) -par1(1)]), lineColor(k), 'LineWidth',2);
            fig = [fig  p1 p2];
        end  
        idx = allArray.NSub == j; 
        par2 = fitDOG(allArray.dDur(idx), allArray.repError(idx));
        plot(xlist, gaussmf(xlist,[par2(2) par2(1)])-gaussmf(xlist,...
                [par2(2) -par2(1)]), lineColor(3), 'LineWidth',2);
        
        plot([-2, 2],[0,0],'k-'); % line
        plot([0, 0],[-2,2],'k-'); % line
        xlabel('Previous trial- current trial (targert duration: secs)');
        ylabel('Error on current trial (secs)');
        %legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'Location','southeast');
  
    end
    hold off;
    set(gca,'fontsize',8)
    sgtitle(strcat('Sequentual effect (RawDat)',{' '}, filename));
    saveas(gcf, strcat('../figures/sequential_raw_NSub_', filename,'.png'));

    
    
     %% grand average     
    figure; 
    fig = [];
    hold on;
    for k = 1:2     
        idx = allArray.target == k; 
        x = allArray.dDur(idx); 
        y = allArray.repError(idx); 
        if (k == 1)
             p1 = plot(x, y, 'r.'); 
        else 
             p1 = plot(x, y, 'g.'); 
        end
      
        par1 = fitDOG(x, y);
        p2 = plot(xlist, gaussmf(xlist,[par1(2) par1(1)])-gaussmf(xlist,...
            [par1(2) -par1(1)]), lineColor(k), 'LineWidth',2);
        fig = [fig  p1 p2];
    end  
    
    %red and green
    par3 = fitDOG(allArray.dDur, allArray.repError);
    p3 = plot(xlist, gaussmf(xlist,[par3(2) par3(1)])-gaussmf(xlist,...
                [par3(2) -par3(1)]), lineColor(3), 'LineWidth',2);
    fig = [fig  p3];
    par_seq = par3;

    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('previous trial -current trial (target duration:secs)');
    ylabel('Error on current trial (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
    title(strcat('Sequentual effect (RawDat)',{' '}, filename));
    saveas(gcf,strcat('../figures/sequential_raw', filename,'.png'));
     
end 

%% visualization to show sequential effect dDur = TD_(2,n-1) - TD_(2,n)
function par_seq = plotSequentialEffect(allArray, filename)
    lineColor = ['r', 'g', 'b'];
    dDur_dat = allArray(allArray.dDurRound~=0,:);
    avg_dDur = grpstats(dDur_dat, {'NSub','dDurRound','target'},{'mean',...
        'sem'}, 'DataVars', {'repError','RP'});

    nsub_dDur = length(unique(avg_dDur.NSub));
    xlist = min(avg_dDur.dDurRound)-2: 0.05 :max(avg_dDur.dDurRound)+2;
    figure;

    for j = 1:nsub_dDur 
        fig = [];
        subplot(4,4,j);
        hold on;
        for k = 1:2     
            idx = avg_dDur.NSub == j & avg_dDur.target == k; 
            x = avg_dDur.dDurRound(idx); 
            y = avg_dDur.mean_repError(idx); 
           
            if (k == 1)
                 p1 = plot(x, y, 'r.'); 
            else 
                 p1 = plot(x, y, 'g.'); 
            end
            
            par1 = fitDOG(x, y);
            p2 = plot(xlist, gaussmf(xlist,[par1(2) par1(1)])-gaussmf(xlist,...
                [par1(2) -par1(1)]), lineColor(k), 'LineWidth',2);
            fig = [fig  p1 p2];
        end  
        idx = avg_dDur.NSub == j; 
        par2 = fitDOG(avg_dDur.dDurRound(idx), avg_dDur.mean_repError(idx));
        plot(xlist, gaussmf(xlist,[par2(2) par2(1)])-gaussmf(xlist,...
                [par2(2) -par2(1)]), lineColor(3), 'LineWidth',2);
        
        plot([-2, 2],[0,0],'k-'); % line
        plot([0, 0],[-2,2],'k-'); % line
        xlabel('previous trial -current trial (target duration:secs)');
        ylabel('Error on current trial (secs)');
        %legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'Location','southeast');
  
    end
    hold off;
    set(gca,'fontsize',8)
    sgtitle(strcat('Sequentual effect',{' '}, filename));
    saveas(gcf, strcat('../figures/sequential_NSub', filename,'.png'));

    
    
     %% grand average     
    all_avg_dDur = grpstats(avg_dDur, {'target','dDurRound'},{'mean','sem'}, ...
        'DataVars', {'mean_repError','mean_RP'});
    figure; 
    fig = [];
    hold on;
    for k = 1:2     
        idx = all_avg_dDur.target == k; 
        x = all_avg_dDur.dDurRound(idx); 
        y = all_avg_dDur.mean_mean_repError(idx); 
        if (k == 1)
             p1 = plot(x, y, 'r.'); 
        else 
             p1 = plot(x, y, 'g.'); 
        end
      
        par1 = fitDOG(x, y);
        p2 = plot(xlist, gaussmf(xlist,[par1(2) par1(1)])-gaussmf(xlist,...
            [par1(2) -par1(1)]), lineColor(k), 'LineWidth',2);
        fig = [fig  p1 p2];
    end  
    
    %red and green
    par3 = fitDOG(all_avg_dDur.dDurRound, all_avg_dDur.mean_mean_repError);
    p3 = plot(xlist, gaussmf(xlist,[par3(2) par3(1)])-gaussmf(xlist,...
                [par3(2) -par3(1)]), lineColor(3), 'LineWidth',2);
    fig = [fig  p3];
    par_seq = par3;

    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('previous trial-current trial (target duration:secs)');
    ylabel('Error on current trial (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
    title(strcat('Sequentual effect',{' '}, filename));
    saveas(gcf,strcat('../figures/sequential', filename,'.png'));
     
end 



 %% visualization of raw data to show contrast effect 
 %NTD_n -  TD_n (TD:target duration; NTD: not target duration) 
function par_con = ShowContrastEffectRawDat(allArray, filename)
    nsub = length(unique(allArray.NSub));
    lineColor = ['r', 'g', 'b'];
    xlist = min(allArray.nttDur): 0.05 :max(allArray.nttDur);
    
    figure;
    
    for j = 1:nsub 
        subplot(4,4,j);
        hold on;      
        for k = 1:2     
            idx = allArray.NSub == j & allArray.target == k; 
            x = allArray.nttDur(idx); 
            y = allArray.repError(idx); 
            if (k == 1)
                 p1 = plot(x, y, 'r.'); 
            else 
                 p1 = plot(x, y, 'g.'); 
            end
            
            par = fitDOG(x, y);
            p2 = plot(xlist, gaussmf(xlist,[par(2) par(1)])-gaussmf(xlist,...
                [par(2) -par(1)]), lineColor(k), 'LineWidth',2);
          
        end     
        idx = allArray.NSub == j;
        par2 = fitDOG(allArray.nttDur(idx), allArray.repError(idx));
        p3 = plot(xlist, gaussmf(xlist,[par2(2) par2(1)])-gaussmf(xlist,...
                [par2(2) -par2(1)]), lineColor(3), 'LineWidth',2);
        plot([-2, 2],[0,0],'k-'); % line
        plot([0, 0],[-2,2],'k-'); % line
        xlabel('nontarget - target duration (secs)');
        ylabel('Error on current trial (secs)');
        legend('Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
   
    end
    hold off;
    set(gca,'fontsize',8)
    sgtitle(strcat('Contrast Effect(Raw data)',{' '}, filename));
    saveas(gcf,strcat('../figures/ContrastEffect_raw_NSub', filename,'.png'));
    
     %% grand average  
    figure; 
    fig = [];
    hold on;
    for k = 1:2     
        idx =  allArray.target == k; 
        x = allArray.nttDur(idx); 
        y = allArray.repError(idx); 
        if (k == 1)
             p1 = plot(x, y, 'r.'); 
        else 
             p1 = plot(x, y, 'g.'); 
        end
        
        par = fitDOG(x, y);
        p2 = plot(xlist, gaussmf(xlist,[par(2) par(1)])-gaussmf(xlist,...
        [par(2) -par(1)]), lineColor(k), 'LineWidth',2);
        fig = [fig  p1 p2];
    end    
     %red and green
    par3 = fitDOG(allArray.nttDur, allArray.repError);
    plot(xlist, gaussmf(xlist,[par3(2) par3(1)])-gaussmf(xlist,...
                [par3(2) -par3(1)]), lineColor(3), 'LineWidth',2);
    fig = [fig  p3];
    par_con = par3;
    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('nontarget - target (secs)');
    ylabel('Error on current trial (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
    title(strcat('Contrast Effect(RawDat)',{' '}, filename));
    saveas(gcf,strcat('../figures/ContrastEffect', filename,'.png'));
    
end


 %% visualization of rund date to show contrast effect 
     % TD_n - NTD_n (TD:target duration; NTD: not target duration) 
function par_con = ShowContrastEffect(allArray, filename)

    avg_nttDur = grpstats(allArray, {'NSub','nttDur_Round','target'},{'mean','sem'},...
        'DataVars', {'repError','RP'});
    all_avg_nttDur = grpstats(avg_nttDur, {'target','nttDur_Round'},{'mean','sem'}, ...
        'DataVars', {'mean_repError','mean_RP'});
    nsub_nttDur = length(unique(avg_nttDur.NSub));
    lineColor = ['r', 'g', 'b'];
    xlist = min(avg_nttDur.nttDur_Round)-2: 0.05 :max(avg_nttDur.nttDur_Round)+2;
    
    figure;
    
    for j = 1:nsub_nttDur 
        subplot(4,4,j);
        hold on;      
        for k = 1:2     
            idx = avg_nttDur.NSub == j & avg_nttDur.target == k; 
            x = avg_nttDur.nttDur_Round(idx); 
            y = avg_nttDur.mean_repError(idx); 
            if (k == 1)
                 p1 = plot(x, y, 'r.'); 
            else 
                 p1 = plot(x, y, 'g.'); 
            end
            
            par = fitDOG(x, y);
            p2 = plot(xlist, gaussmf(xlist,[par(2) par(1)])-gaussmf(xlist,...
                [par(2) -par(1)]), lineColor(k), 'LineWidth',2);
          
        end     
        idx = avg_nttDur.NSub == j;
        par2 = fitDOG(avg_nttDur.nttDur_Round(idx), avg_nttDur.mean_repError(idx));
        p3 = plot(xlist, gaussmf(xlist,[par2(2) par2(1)])-gaussmf(xlist,...
                [par2(2) -par2(1)]), lineColor(3), 'LineWidth',2);
        plot([-2, 2],[0,0],'k-'); % line
        plot([0, 0],[-2,2],'k-'); % line
        xlabel('nontarget - target duration (secs)');
        ylabel('Error on current trial (secs)');
        legend('Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
   
    end
    hold off;
    set(gca,'fontsize',8)
    sgtitle(strcat('Contrast Effect',{' '}, filename));
    saveas(gcf,strcat('../figures/ContrastEffect_NSub', filename,'.png'));
    
     %% grand average  
    figure; 
    fig = [];
    hold on;
    for k = 1:2     
        idx =  all_avg_nttDur.target == k; 
        x = all_avg_nttDur.nttDur_Round(idx); 
        y = all_avg_nttDur.mean_mean_repError(idx); 
        if (k == 1)
             p1 = plot(x, y, 'r.'); 
        else 
             p1 = plot(x, y, 'g.'); 
        end
        
        par = fitDOG(x, y);
        p2 = plot(xlist, gaussmf(xlist,[par(2) par(1)])-gaussmf(xlist,...
        [par(2) -par(1)]), lineColor(k), 'LineWidth',2);
        fig = [fig  p1 p2];
    end    
     %red and green
    par3 = fitDOG(all_avg_nttDur.nttDur_Round, all_avg_nttDur.mean_mean_repError);
    plot(xlist, gaussmf(xlist,[par3(2) par3(1)])-gaussmf(xlist,...
                [par3(2) -par3(1)]), lineColor(3), 'LineWidth',2);
    fig = [fig  p3];
    par_con = par3;
    plot([-2, 2],[0,0],'k-'); % line
    plot([0, 0],[-2,2],'k-'); % line
    hold off;
    xlabel('nontarget - target (secs)');
    ylabel('Error on current trial (secs)');
    legend(fig, 'Red target','Red(fit)', 'Green target', 'Green(fit)', 'All (fit)','Location','northeast');
    title(strcat('Contrast Effect',{' '}, filename));
    saveas(gcf,strcat('../figures/ContrastEffect', filename,'.png'));
    
end

function par = fitDOG(x, y)
     ft = fittype('1/(sig*sqrt(2*pi)) * exp(-(x-mu)^2/(2*sig^2))- 1/(sig*sqrt(2*pi)) * exp(-(x+mu)^2/(2*sig^2))',...
        'independent', 'x','coefficients', {'mu','sig'});
     [obj,gof,opt]  = fit(x, y,ft)
     par = coeffvalues(obj); 

end