function y=datainput_direct(xx,prob) %objective function
x=xx'; % make sure vector is row vector

%test
% x=[2 3 4 5 2 3 4 5 1 1 2 3 3 1 2 3 4 5 1 1 1 1 1 1 1 1 1 1 1 1 1];
% x=[2:20 1 1 2 3 3 1 1 2 3 3 1 1 2 3 3 1 1 2 3 3 2 3 3 ];
%test


n=prob.n;
n_S=prob.n_S;
n_M=prob.n_M;
A=zeros(n,n_S);
y=0;
C=zeros(n,n_S);
for j=1:n_S
    job_seq=[0,1,0];
    for i=1:n-1
        job_seq=[job_seq(1:x((j-1)*(n-1)+i)) i+1 job_seq(x((j-1)*(n-1)+i)+1:end)];
    end
    job_seq=job_seq(2:n+1);
    assign=x(n_S*(n-1)+(j-1)*n+1:n_S*(n-1)+j*n);
    pr_speed=x((n-1)*n_S+n*n_S+sum(n_M(1:j-1))+1:(n-1)*n_S+n*n_S+sum(n_M(1:j)));
    ST_rate=x(prob.int_dim+(j-1)*n+1:prob.int_dim+j*n);
    [stage_C,stage_EP,stage_ETI]=timingByRate(job_seq,assign,pr_speed,A(:,j),j,ST_rate,prob);
    y=y+stage_EP+stage_ETI;
    C(:,j)=stage_C;

    
    if j<n_S
        l=prob.l;
        E_Tr=zeros(n,1);
        tr_speed=x((n-1)*n_S+n*n_S+sum(n_M)+(j-1)*n+1:(n-1)*n_S+n*n_S+sum(n_M)+j*n);
        assign_next=x(n_S*(n-1)+j*n+1:n_S*(n-1)+(j+1)*n);
        for i=1:n
            A(i,j+1)=C(i,j)+l(j,assign(i),assign_next(i),tr_speed(i));
            E_Tr(i)=l(j,assign(i),assign_next(i),tr_speed(i))*prob.e_t(tr_speed(i));
        end
        y=y+sum(E_Tr);
    end
    
end







end %
