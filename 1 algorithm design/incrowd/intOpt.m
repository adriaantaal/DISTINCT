function X = intOpt(Asmaller,Y,Lambda,xo,interiorSolveOption)
%  Gives an exact basis pursuit answer to the dense problem handed to it.
%  Lambda,
switch interiorSolveOption
    case 1 % Use quadprog
        H = Asmaller' * Asmaller;
        LambdaVect = 2 * Lambda * ones(size(Asmaller,2),1);
        f = ((LambdaVect./2)-Asmaller' * (Y));

        if size(Asmaller,2) > size(Asmaller,1)
%             disp('Using qps_mq');
            %  There is a bug in MATLAB's native quadprog, such that it
            %  does not respect solution constraints if size(Asmaller,2) > size(Asmaller,1)

            %  This uses SPM's quadprog routine instead:  see
            %  http://sigpromu.org/quadprog for details
            [X,err] = qps_mq(H,f,zeros(size(Asmaller,2),1),inf(size(Asmaller,2),1));
            if err > 0
                error('Inexact solution found.');
            end
        else
            A = -1*eye(size(Asmaller,2));%; eye(size(Asmaller,2))];
            b = zeros(1,length(f));% maxVal*ones(1,length(f))];
            %OPTIONS = optimget
            %OPTIONS = optimset('Display','iter','LargeScale','off');
            maxiter = 1e5;
            [X,fval,exitflag,output] = quadprog(H,f,A,b,[],[],[],[],xo,optimset('MaxIter',maxiter,'Display','off','LargeScale','off'));
%                 disp(['exitflag is ' num2str(exitflag)])    
                              
            if exitflag < 0
                disp(['exitflag is ' num2str(exitflag)])
%                 maxiter = maxiter/10;
%                 disp(['running ' num2str(maxiter) ' max iterations'])
%                 disp(['Trying again ' num2str(maxiter) ' max iterations'])
%                 [X,fval,exitflag] = quadprog(H,f,A,b,[],[],[],[],xo,optimset('MaxIter',maxiter,'Display','off'));

                disp('trying l1_ls_nonneg')
                
                [X,status,~] = l1_ls_nonneg(Asmaller,Y,Lambda,1e-3,1);
                
                if strcmp(status,'Failed')
                    disp(['trying SolveBP'])
                    X = SolveBP(Asmaller, Y, size(Asmaller,2), inf, Lambda);
                end
            end
        end
        
    case 2 % Use homotopy
        maxiter = 1e5;
        X = BPDN_positivity_function(Asmaller,Y,Lambda, maxiter);
    case 3
        X = l1_ls_nonneg(Asmaller,Y,Lambda,1e-3,1);
    case 4
        X = SolveBP(Asmaller, Y, size(Asmaller,2), inf, Lambda);
end
