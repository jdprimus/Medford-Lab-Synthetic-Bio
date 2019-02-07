function status = interruptFun(t,y,flag,interrupt_time)   %#ok<INUSL>
% function to interrupt ODE solver if integration takes too long
    persistent INIT_TIME;
    status = 0;
    switch(flag)
        case 'init'
            INIT_TIME = tic;
        case 'done'
            clear INIT_TIME;
        otherwise
            elapsed_time = toc(INIT_TIME);
            if elapsed_time > interrupt_time
                status = 1;
                clear INIT_TIME;
                str = sprintf('%.6f',elapsed_time);
                error('interruptFun:Interrupt',...
                     ['Interrupted integration. Elapsed time is ' str ' seconds.']);
            end
%             if NonNeg == 1
%                 if sum(y(2:5) < 0.1) > 0
%                 status = 1;
%                 str = sprintf('%.6f',elapsed_time);
%                 format long
%                 disp(y)
%                 error('interruptFun:Interrupt',...
%                      ['Interrupted integration. Values less than zero recorded. ' str ' seconds.' ]);
%                 end
%             end
    end
end