M=100;                 
D=1;                  
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];

post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
s=[1*ones(Ne,M);-1*ones(Ni,M)];
sd=zeros(N,M);                  

delays = cell(N,D);
for i=1:Ne
    p=randperm(N);
    post(i,:)=p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  
    end;
end;
for i=Ne+1:N
    p=randperm(Ne);
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    
end;
post(111,:)=randperm(Ne,M);
post(222,:)=randperm(Ne,M);

pre = cell(N,1);
aux = cell(N,1);
for i=1:Ne
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end
    end
end
  
STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings
%selecting excitatory neurons for action representation
n1_right=post(222,:);
n1_left=post(111,:);
S1_neurons=111;
S2_neurons=222;

save('intial.mat','pre','delays','post','aux');

%4 inputs representation
I_input=12*(ones(N,4));
for i=1:1000
  if ~any(i==S1_neurons)
    I_input(i,1)=0; %earlier kept 0
  end
  if ~any(i==S2_neurons)
    I_input(i,2)=0;
  end
   if ~any(i==333)
    I_input(i,3)=0;
   end
 
   I_input(i,4)=0; %S4 =0 
end

T=5000;        
DA=0;         
rew=[];
neg_rew=[];
stimulus_time=0;
alt=0;
input=1;
count_left=0;
count_right=0;
count_no_action=0;
incorrect_L=0;
incorrect_R=0;
window=100;
record_mat=[];

for sec=1:T  
     fired_R=0;
     fired_L=0;
     fired_S1=0;
     fired_S2=0;
     fires=zeros(N,1);
     avg_firing=0; 
    if mod(sec,5000)==1
        save(("weight_noise_"+ sec +".mat"),'s');
    end
  for t=1:1000                         
        if(mod(sec,5)==0)
            if(t==1)
                sd=zeros(N,M);                     
                DA=0;
                input=randperm(2,1);
                input2= randperm(2,1)+2;
                I=I_input(:,input)+I_input(:,input2);
                stimulus_time=sec*1000+t;
                %fprintf('sec:%d\t correct L:%d \tR:%d \tNA:%d \n%d\t Incorrect L:%d R:%d dw>1:%d\n',sec,count_left,count_right,count_no_action,input,incorrect_L,incorrect_R,dw_greater_than_one);
            else
                if(t>=window)
                   I=zeros(N,1);
                   I(ceil(rand*N))=10;
                end
            end
        else
            I=zeros(N,1);
            I(ceil(rand*N))=10;
        end
        
    fired = find(v>=30);                
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
      k=k-1;
    end
    
    v=v+0.5*((0.04*v+5).*v+140-u+ I);   
   % v=v+0.5*((0.04*v+5).*v+140-u+ I);    
    u=u+a.*(0.2*v-u);                  
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     
    
    if(t<=window+50 && mod(sec,5)==0)
        temp_R=0;
        temp_L=0;
        temp_s1=0;
        temp_s2=0;
        temp_fires=zeros(N,1);
        avg_firing=avg_firing+size(fired,1);
        for i=1:size(fired,1)
            if any(fired(i)==n1_right)
                temp_R=temp_R+1;
            end
            if any(fired(i)==n1_left)
                temp_L=temp_L+1;
            end
            if any(fired(i)==S1_neurons)
                temp_s1=temp_s1+1;
            else
                if any(fired(i)==S2_neurons)
                    temp_s2=temp_s2+1;
                end
            end
            temp_fires(fired(i))=temp_fires(fired(i))+1;
        end
        
      
        fired_S1=fired_S1+temp_s1;
        fired_S2=fired_S2+temp_s2;
        fired_R=fired_R+temp_R;
        fired_L=fired_L+temp_L;
        fires=fires+temp_fires;
        if (t==window+50)  
            action_executed=[0,0,0]; %[left,right,no action]
            correct=0; %correct=1 for correct case and correct=0 for wrong case  
            avg_firing=avg_firing;
            %fprintf('size: avg: %d fired_ L: %d, R: %d \n',avg_firing,fired_L,fired_R);
            reward_time=-999;
                             
            if(fired_L>fired_R)
                action_executed=[1,0,0];
                reward_time=sec*1000+t+floor(1000/(fired_L-fired_R));      
                if (input==1) 
                    count_left=count_left+1;
                    correct=1;
                    %fprintf('action:left\tcount:%d \n',count_left);
                    %reward
                    rew=[rew,reward_time];   %delay of 1 second
                else
                    incorrect_L=incorrect_L+1;
                    neg_rew=[neg_rew,reward_time];
                end
    
            elseif(fired_R>fired_L)
                action_executed=[0,1,0];
                reward_time=sec*1000+t+floor(1000/(fired_R-fired_L));           
                
                if (input==2)
                  correct=1;
                  count_right=count_right+1;
                    %fprintf('action:right \t  count:%d \n',count_right);
                    %reward
                  
                    rew=[rew,reward_time];   %delay of 1 second
                else
                     incorrect_R=incorrect_R+1;
                     neg_rew=[neg_rew,reward_time];
                end
            else 
                 action_executed=[0,0,1];
                 count_no_action=count_no_action+1;
                %fprintf('action:NA \t count:%d \n',count_no_action);
            end
            record=[correct,input,input2,action_executed,avg_firing,fired_L,fired_R,fired_S1,fired_S2];
            for k=1:N
                record=[record,fires(k)];
            end
            record_mat=[record_mat;record];
        end
    end
  
  if any(rew==sec*1000+t)
      DA=0.5;
  end
   
%   if any(neg_rew==sec*1000+t)
%       DA=-0.5;
%   end  
  
  DA=DA*0.99;
    if (mod(t,10)==0)
        s(1:Ne,:)=max(0,s(1:Ne,:)+(DA)*sd(1:Ne,:));
        sd=0.99*sd;
    end
  
  end

%   subplot(2,1,1)
%   plot(firings(:,1),firings(:,2),'.');
%   axis([0 1000 0 N]); 
% %   subplot(2,2,3);
% %   plot(0.001*(1:(sec*1000+t)), shist(1:sec*1000+t,2),0.001*rew,0*rew,'rx');
%   subplot(2,2,4);
%   hist(s(find(s>0)),60*(0.01:0.01:1)); % only excitatory synapses
%   drawnow;

  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
 
end
save("data.mat",'record_mat');
save("value_noise_.mat",'STDP','s','sd');
