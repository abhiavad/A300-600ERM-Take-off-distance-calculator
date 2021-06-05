function [s_ga_m,s_rot_m,s_climb_m,s,Tb,v_r_m_s]=takeoffdist_m(alt_km,wto_ratio)
g=9.81;
W_to_kn=1685;
W_to_n=W_to_kn*1000;
w_to_given_n=W_to_n*wto_ratio;

alt_m=alt_km*1000;
[~,~,~,rho_kg_m3]=atmosisa(alt_m);
[~,~,~,rho0_kg_m3]=atmosisa(0);
sigma=rho_kg_m3/rho0_kg_m3;

T_kn=275*2*sigma;
T_n=T_kn*1000;

S_m2=260;
AR=7.73;
e=0.91;

cd_0=0.012;
cl_max=1.1;
cl_climb=0.9*cl_max;
delta_cd_0_lg=0.004;
delta_cd_0_flap=0.003;
cl_gr=0.1;
mu=0.02;

k=1/pi*e*AR;
t_r_s=3;

delta_cd_0=delta_cd_0_lg+delta_cd_0_flap;
cd_gr=cd_0+ delta_cd_0+k*cl_gr^2;

v_stall_m_s=sqrt(2*w_to_given_n/(rho_kg_m3*S_m2*cl_max));
v_r_m_s=1.2*v_stall_m_s;


i=1;
v_i=zeros(500,1);
q_i=zeros(500,1);
q_i_s=zeros(500,1);
d_i=zeros(500,1);
l_i=zeros(500,1);
a_i=zeros(500,1);
delta_v_i=zeros(500,1);
delta_s_i=zeros(500,1);
v_i_1=zeros(500,1);
s_i_1=zeros(500,1);
t_i_1=zeros(500,1);

s_i=zeros(500,1);
t_i=zeros(500,1);
delta_t=0.25;

while v_i(i)<v_r_m_s

q_i(i)=0.5*rho_kg_m3*(v_i(i))^2;
q_i_s(i)=q_i(i)*S_m2;
d_i(i)=cd_gr*q_i_s(i);
l_i(i)=cl_gr*q_i_s(i);
a_i(i)=(T_n-d_i(i)-mu*(w_to_given_n-l_i(i)))*g/w_to_given_n;
delta_v_i(i)=a_i(i)*delta_t;
delta_s_i(i)=v_i(i)*delta_t+0.5*a_i(i)*(delta_t^2);
v_i(i+1)=v_i(i)+delta_v_i(i);
v_i_1(i)=v_i(i+1);
s_i(i+1)=s_i(i)+delta_s_i(i);
s_i_1(i)=s_i(i+1);
t_i(i+1)=t_i(i)+delta_t;
t_i_1(i)=t_i(i+1);

i=i+1;
end

%Tb=table(v_i(1:i),q_i_s(1:i),d_i(1:i),l_i(1:i),a_i(1:i),delta_v_i(1:i),delta_s_i(1:i),v_i_1(1:i),s_i_1(1:i),t_i_1(1:i));
Tb=table(v_i,q_i_s,d_i,l_i,a_i,delta_v_i,delta_s_i,v_i_1,s_i_1,t_i_1);
Tb=Tb(1:i-1,:);
Tb.Properties.VariableNames={'Velocity (Vi)','Dynamic Head (qiS)','Drag (Di)','Lift (Li)','Acceleration (ai)','Change in Velocity (ΔVi)','Delta Position(Δsi)','Velocity at (i+1)Δt(V(i+1))','Position (i+1)Δt((s(i+1))','Next time step (t(i+1))'};
s_ga_m=s_i(i);
s_rot_m=3*v_r_m_s;
sin_g=(T_n-((cd_0+delta_cd_0+k*cl_climb^2)*(0.5*rho_kg_m3*v_r_m_s^2)))/w_to_given_n;
gamma=asin(sin_g);
tan_g=tan(gamma);
s_climb_m=10.7/tan_g;
s=s_ga_m+s_rot_m+s_climb_m;

end