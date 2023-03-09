clear all
comnum=50;    %number of components

%number of samples in weighted random sampling
MCnum1=10000;


%definition of component failure probability
alpha1=unifrnd(0.8,1.6,3,comnum);
beta1=unifrnd(1.3,2.1,3,comnum);

%alpha1=unifrnd(1.2,1.8,3,comnum);
%beta1=unifrnd(2.3,2.9,3,comnum);


rdc=0;

%1¡¢calculation of minimum divided unit
s1=comnum;
cellnum=0;
judgem1=0;
while s1>4
    if mod(s1,2) == 0
        s1=s1/2;
    else
        s1=(s1-1)/2;
        judgem1=1;
    end
    cellnum=cellnum+1;
end
if judgem1==1
    s2=s1+1;
else
    s2=s1;
end
n2=comnum-s1*2^cellnum;
n1=2^cellnum-n2;


SS=[];

for lll=1:25
    fnum=lll;
    hhh=1;
    for T=0.1:0.1:2
        rdc=0;
        for n=1:comnum
            reliability(n)=1-gamcdf(T,alpha1(1,n),beta1(1,n));
            %reliability(n)=exp(-(T/alpha1(1,n))^beta1(1,n));
        end
        
        for s=cellnum:-1:0
            subsystemN=[];
            if s==cellnum
                for t=1:2^s
                    subsystemR=1;
                    if t<=n1
                        y=1:s1;
                        subsystemN=zeros(1,s1+1);
                        for m=(s1*(t-1)+1):(s1*t)
                            subsystemR=reliability(m)*subsystemR;
                        end
                        for n=1:fnum
                            C13=nchoosek(y,n);
                            subsystem=zeros(1,length(C13(:,1)));
                            for k=1:length(C13(:,1))
                                subsystemR1=subsystemR;
                                for h=1:n
                                    subsystemR1=subsystemR1*(1-reliability(C13(k,h)+s1*(t-1)))/reliability(C13(k,h)+s1*(t-1));
                                end
                                subsystem(k)=subsystemR1;
                            end
                            subsystemN(n+1)=sum(subsystem);
                        end
                        subsystemN(1)=subsystemR;
                        A{s+1,t}=subsystemN;
                    else
                        y=1:s2;
                        subsystemN=zeros(1,s2+1);
                        for m=(n1*s1+s2*(t-n1-1)+1):(n1*s1+s2*(t-n1))
                            subsystemR=reliability(m)*subsystemR;
                        end
                        for n=1:fnum
                            C12=nchoosek(y,n);
                            subsystem=zeros(1,length(C12(:,1)));
                            for k=1:length(C12(:,1))
                                subsystemR1=subsystemR;
                                for h=1:n
                                    subsystemR1=subsystemR1*(1-reliability(C12(k,h)+n1*s1+s2*(t-n1-1)))/reliability(C12(k,h)+n1*s1+s2*(t-n1-1));
                                end
                                subsystem(k)=subsystemR1;
                            end
                            subsystemN(n+1)=sum(subsystem);
                        end
                        subsystemN(1)=subsystemR;
                        A{s+1,t}=subsystemN;
                    end
                end
            else
                for g=1:2^s
                    subsystemg=zeros(1,fnum);
                    lastlayer1=A{s+2,(2*g-1)};
                    lastlayer2=A{s+2,2*g};
                    for f=0:fnum
                        layerR=0;
                        if f<=length(lastlayer2)-1
                            for u=0:f
                                layerR=layerR+lastlayer1(u+1)*lastlayer2(f-u+1);
                            end
                        else
                            for u=(f-length(lastlayer2)+1):length(lastlayer2)-1
                                layerR=layerR+lastlayer1(u+1)*lastlayer2(f-u+1);
                            end
                        end
                        subsystemg(f+1)=layerR;
                    end
                    A{s+1,g}=subsystemg;
                end
            end
        end
        
        failstate=0;
        for n=1:MCnum1
            state=ones(1,comnum);
            B=[];
            B(1,1)=fnum;
            for i=1:cellnum
                for j=1:2^(i-1)
                    goalfnum=B(i,j);
                    generate1=A{i+1,2*j-1};
                    generate2=A{i+1,2*j};
                    w=[];
                    for k=1:goalfnum+1
                        w(k)=generate1(k)*generate2(goalfnum+2-k);
                    end
                    random=rand(1,length(w));
                    maxkey=0;
                    for k=1:goalfnum+1
                        if w(k)~=0
                            key(k)=random(k)^(sum(w)/w(k));
                            if key(k)>maxkey
                                maxkey=key(k);
                                maxindex=k-1;
                            end
                        end
                    end
                    B(i+1,2*j-1)=maxindex;
                    B(i+1,2*j)=goalfnum-maxindex;
                end
            end
            %µ×²ã³éÑù
            for i=1:n1
                if B(cellnum+1,i)~=0
                    rall=1;
                    for j=1:s1
                        rall=reliability(s1*(i-1)+j)*rall;
                    end
                    generatemode=nchoosek(1:s1,B(cellnum+1,i));
                    probability=[];
                    for k=1:length(generatemode(:,1))
                        for s=1:length(generatemode(1,:))
                            probability(k)=rall/reliability(s1*(i-1)+generatemode(k,s))*(1-reliability(s1*(i-1)+generatemode(k,s)));
                        end
                    end
                    
                    random2=rand(1,length(generatemode(:,1)));
                    maxkey2=0;
                    for k=1:length(generatemode(:,1))
                        key2(k)=random2(k)^(sum(probability)/probability(k));
                        if key2(k)>maxkey2
                            maxkey2=key2(k);
                            maxindex2=k;
                        end
                    end
                    for s=1:length(generatemode(1,:))
                        state(generatemode(maxindex2,s)+s1*(i-1))=0;
                    end
                end
            end
            
            for i=1:n2
                if B(cellnum+1,i+n1)~=0
                    rall=1;
                    for j=1:s2
                        rall=reliability(s1*n1+s2*(i-1)+j)*rall;
                    end
                    generatemode=nchoosek(1:s2,B(cellnum+1,i+n1));
                    probability=[];
                    for k=1:length(generatemode(:,1))
                        for s=1:length(generatemode(1,:))
                            probability(k)=rall/reliability(s1*n1+s2*(i-1)+generatemode(k,s))*(1-reliability(s1*n1+s2*(i-1)+generatemode(k,s)));
                        end
                    end
                    
                    random2=rand(1,length(generatemode(:,1)));
                    maxkey2=0;
                    for k=1:length(generatemode(:,1))
                        key2(k)=random2(k)^(sum(probability)/probability(k));
                        if key2(k)>maxkey2
                            maxkey2=key2(k);
                            maxindex2=k;
                        end
                    end
                    for s=1:length(generatemode(1,:))
                        state(generatemode(maxindex2,s)+s1*n1+s2*(i-1))=0;
                    end
                end
            end
            
            %ÅÐ¶ÏÃ¿¸ö×´Ì¬º¯Êý
            for js=1:2:comnum-1
                if state(js)==0 && state(js+1)==0
                    failstate=failstate+1;
                    break
                end
            end
        end
        SS(lll,hhh)=1-failstate/MCnum1;
        
        hhh=hhh+1;

    end
end

%analytical solution
fiidd=[];
hhh=1;
for lll=1:5
    fnum=test(lll);
    zz1=1;
    zz2=1;
    if fnum~=0
        for pp=1:fnum
            zz1=zz1*(50-pp+1);
            zz2=zz2*(50-2*pp+2);
        end
    end
    fiidd(hhh)=zz2/zz1;
    
    hhh=hhh+1;
end