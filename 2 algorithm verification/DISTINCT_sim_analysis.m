%DISTINCT_sim_analysis
%

names = {'DISTINCT' , ...
         'ICOpt' , ...
         'SL0' , ...
         'NIHT' , ...
         'DSGD' , ...
         'DINV'}

%% do fitting on the time
clear N D yy p

figure; plot(L,tcomp);
set(gca,'xscale','log');

figure; set(gcf,'Color',[1 1 1]);
for jj = 1:size(tcomp,2)
%     p(:,jj) = polyfit(sqrt(L),tcomp(:,jj),1)
%     yy(:,jj) = polyval(p(:,jj),sqrt(L));

    p(:,jj) = polyfit(log(L),tcomp(:,jj),1)
    yy(:,jj) = polyval(p(:,jj),log(L));

%     [N(:,jj),D(:,jj)] = ratpolyfit(L,tcomp(:,jj),2,2);
%     yy(:,jj)=polyval(N(:,jj),L)./polyval(D(:,jj),L);

%     myfittype=fittype('a +b*log(x)','dependent', {'y'}, 'independent',{'x'},'coefficients',{'a,b'});
%     myfit=fit(L',tcomp(:,jj)',myfittype,'StartPoint',[1 1]);

    subplot(2,3,jj); hold on;
%     plot(L,tcomp(:,jj));
%     plot(L,yy(:,jj));
    loglog(L,tcomp(:,jj));
    loglog(L,yy(:,jj));
    title(names{jj});
    set(gca,'xscale','log');
    set(gca,'yscale','log');
end


