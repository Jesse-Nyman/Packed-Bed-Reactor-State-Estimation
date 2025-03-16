function  Parsave(j,n,xout,yout,Ut,theta)
    save(sprintf('%s%d_%d_%s','..\Data\',j,n,'xout.txt'),'xout', '-ascii');
    save(sprintf('%s%d_%d_%s','..\Data\',j,n,'yout.txt'),'yout','-ascii');
    save(sprintf('%s%d_%d_%s','..\Data\',j,n,'Ut.txt'),'Ut','-ascii');
    dlmwrite(sprintf('%s%d_%d_%s','..\Data\',j,n,'theta.txt'),theta);

%     save(sprintf('%s%d_%d_%s','/scratch/work/forsstn1/Data/',j,n,'xout.mat'),'xout');
%     save(sprintf('%s%d_%d_%s','/scratch/work/forsstn1/Data/',j,n,'yout.mat'),'yout');
%     save(sprintf('%s%d_%d_%s','/scratch/work/forsstn1/Data/',j,n,'Ut.mat'),'Ut');
%     save(sprintf('%s%d_%d_%s','/scratch/work/forsstn1/Data/',j,n,'theta.mat'),'theta');
end
