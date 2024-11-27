    %Anil Cengiz - Sep 12,2024
    %This function generates time courses of drug concentration based on the PK
    %parameters of a drug 
    
    function [time,conc]=conc_time_courses_mod(par_struct,pl_ind,tau2,tirz_dose,switch_time)
    tau1=7*24; %hours - Original dosing interval 
    esc_time=switch_time;

    %For Semaglutide 
    if par_struct.name=="semaglutide" 
    ka=par_struct.ka;
    V=par_struct.V;
    ke=par_struct.ke;

    F=par_struct.F; 
    weekly_dose=par_struct.weekly_dose;
    
    %switch_time=80*7*24; %When the dosing regimen is switched to an alternative dosing regimen
    dosing_times_1=0:tau1:(switch_time-tau1); %Dosing times up to switch
    dosing_times_2=(switch_time-tau1+tau2):tau2:switch_time-tau1+tau2+tau2*186;%Dosing times after switch
    dosing_times=[dosing_times_1 dosing_times_2];
    
    max_time=dosing_times(end);
    dt=1;
    time=0:dt:max_time+tau2-dt;
    c=zeros(size(time)) ; %concentration holder
    g=zeros(size(time)); %drug amount holder 
    D=weekly_dose+zeros(size(dosing_times));

    D(1:4)=0.25;
    D(5:8)=0.5;
    D(9:12)=1;
    D(13:16)=1.7;
    D(17:end)=2.4;

    
    if pl_ind==0 %Placebo from day zero
        D=D*0;
    elseif pl_ind==1 %Placebo after switch 
        D((switch_time/(7*24)+1):end)=0;
    end

    %Loop over dosing times 
    for j=1:length(dosing_times)
        i=dosing_times(j)+1;
        %Update drug amount at absorption site
        g(i)=g(i)+D(j)*F;
        y0=[g(i),c(i)];
        
        %Pick correct dosing interval 
        if j<switch_time/(24*7)
            tau=7*24;
        else
            tau=tau2;
        end

        %Flow ODE to find g and c during 1 dosing cycle 
        tspan=0:dt:tau;
        [t,y] = ode45(@(t,y) odefcn_semag(t,y,V,ka), tspan, y0);
        g(i:i+tau)=y(:,1);
        c(i:i+tau)=y(:,2);
    
    end

    %Create a vector that has the timecourse for the whole drug amount (at absorption site) and drug concentration 
    conc=[g(1:(end-1));c(1:(end-1))];
  conc=conc(2,:)'*1e6/par_struct.molar_mass;

    
    
    
    %For Tirzepatide 
    elseif par_struct.name=="tirzepatide"
    ka=par_struct.ka;
    Vc=par_struct.Vc;
    Vp=par_struct.Vp;
    cl=par_struct.cl;
    q=par_struct.q;
    F=par_struct.F;
    

    %esc_time=80*7*24; %Dosing switch time 
    dosing_times_1=0:tau1:(esc_time-tau1); %Dosing times up to switch
    dosing_times_2=(esc_time-tau1+tau2):tau2:esc_time-tau1+tau2+tau2*186; %Dosing times after switch
    dosing_times=[dosing_times_1 dosing_times_2];
    max_time=dosing_times(end);
    dt=1;
    time=0:dt:max_time+tau2-dt;
    
    c=zeros(size(time)) ; %holder for drug concentration in the central compt.
    cp=zeros(size(time)) ; %holder for drug concentration in the peripheral compt.
    g=zeros(size(time)); %holder for drug concentration at absorption site

    %Tirzepatide uses a dose escalation scheme in which, dose is increased
    %by 2.5 mg every 4 weeks (equivalently, every 4 doses)

    D=2.5+zeros(size(dosing_times));
    D(5:end)=5;
    if tirz_dose==5  %if dose is 5mg
   
    elseif tirz_dose==10  %if dose is 10mg    
    D(9:end)=7.5;
    D(13:end)=10;
    
    elseif tirz_dose==15  %if dose is 15mg
           D(9:end)=7.5;
    D(13:end)=10;
    D(17:end)=12.5;
    D(21:end)=15;
    end
    
    if pl_ind==0 %Placebo from day zero
        D=D*0;
    elseif pl_ind==1
            D((esc_time/(7*24)+1):end)=0; %Placebo after switch  
    
    end
    %Loop over dosing times 
    for j=1:length(dosing_times)
        i=dosing_times(j)+1;
    
        %Update drug at absorption site
        g(i)=g(i)+D(j)*F;
        y0=[g(i),c(i),cp(i)];
        
        %Choose correct dosing interval
        if j<esc_time/(7*24)
        tau=7*24;
        else
        tau=tau2;
        end
        %Flow ODE to find g and c and cp during 1 dosing cycle 
        tspan=0:dt:tau;
        [~,y] = ode45(@(t,y) odefcn_tirz(t,y,Vc,Vp,cl,q,ka), tspan, y0);

    
        g(i:i+tau)=y(:,1); %drug amount
        c(i:i+tau)=y(:,2); %central conc
        cp(i:i+tau)=y(:,3); %peri. conc
    end
    
        %Create a vector that has the timecourse for the whole drug amount
        %(at absorption site) and drug concentration at central and
        %peripheral compartments
    conc=[g(1:(end-1));c(1:(end-1));cp(1:(end-1))];
    conc=conc(2,:)'*1e6/par_struct.molar_mass;
    
    end
    
    
    function dydt = odefcn_tirz(~,y,Vc,Vp,cl,q,ka)
    %function dydt = odefcn(t,y,Vc,Vp,cl_time,cl_vec,q,ka)
    %ind=find(t>=cl_time*7*24,1,'last');
    %cl=cl_vec(ind);
    
      dydt = zeros(3,1);
      dydt(1)=-ka*y(1);
      dydt(2) = ka*y(1)/Vc - cl*y(2)/Vc - q/Vc*y(2) + q/Vp*y(3);
      dydt(3)=  q/Vc*y(2) - q/Vp*y(3);
    end
    
    function dydt = odefcn_semag(~,y,V,ka)
      dydt = zeros(2,1);
      dydt(1)=-ka*y(1);
      dydt(2) = ka*y(1)/V - ke*y(2);
    end
    end