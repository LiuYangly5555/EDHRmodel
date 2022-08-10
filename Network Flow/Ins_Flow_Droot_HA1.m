function obj= Ins_Flow_Droot_HA1(SWC1, SWC2, H1, H2, R, l_s, l_d, ln_s, ln_d, SWCr, SWCs, n, Ks, l, D1, D2)

Ra_r1=R.Ra_r1.*10^4./600;
Ra_r1_2=R.Ra_r1_2.*10^4./600;
Ra_r2=R.Ra_r2.*10^4./600;
Ra_s1=R.Ra_s1.*10^4./600;
Ra_s2=R.Ra_s2.*10^4./600;
Ra_b1=R.Ra_b1.*10^4./600;
Ra_b2=R.Ra_b2.*10^4./600;
Rc_r=R.Rc_r.*10^4./600;
Rc_r_2=R.Rc_r2.*10^4./600;
Rc_s=R.Rc_s.*10^4./600;
Rr_s=R.Rr_s.*10^4./600;
Rr_s_d=R.Rr_s.*10^4./600.*2;
Rr_d=R.Rr_d.*10^4./600;

% set leaf water potential to -1 MPa
Hleaf=-10000;

Ks2 = swc2ks(SWC2, SWCr, SWCs, n, Ks, l);
Rs2 = 1./Ks2;
Rrs1_2=Rs2.*ln_d./2./pi./l_s.*2.*24.*6;
Rrs2 = Rs2.*ln_d./2/pi./l_d.*24.*6;

Ks1 = swc2ks(SWC1, SWCr, SWCs, n, Ks, l);
Rs1 = 1./Ks1;
Rrs1 = Rs1.*ln_s./2./pi./l_s.*24.*6;

a=size(Ra_r1);
a=a(1, 2);

for i=1:1:a

    syms q1 q2 q3 p1 p2 p3;
    eqns = [-q1*(2*Rrs1+2*Rr_s+Ra_r1(1, i))+p1*Rc_r(1, i)+q2*(2*Rrs1+2*Rr_s+Ra_r1(1, i))==0,...
        -(q1+p1)*(Ra_r1(1, i)+Ra_s1/2)+p3*Rc_s+(q3-p2)*(Ra_s2/2+Ra_r2(1, i)/2)-p2*Rc_r_2(1, i)+(q2-p1)*(Ra_r1(1, i)+Ra_r1_2(1, i)/2)-p1*Rc_r(1, i)==0,...
        -q3*(Rrs2+Rr_d+Ra_r2(1, i)/2)-p2*Rc_r_2(1, i)-(q2-p1+p2)*(Rrs1_2+Rr_s_d+Ra_r1_2(1, i)/2)==0,...
        (q3-p2-p3)*(Ra_s2/2+Ra_b2)-p3*Rc_s-(q1+p1+p3)*(Ra_s1/2+Ra_b1)==0,...
        (H2-D1-D2/2)-q3*(Rrs2+Rr_d+Ra_r2(1, i)/2)-(q3-p2)*(Ra_r2(1, i)/2+Ra_s2/2)-(q3-p2-p3)*(Ra_s2/2+Ra_b2)-Hleaf==0,...
        (H2-D1-D2/2)+(q2-p1+p2)*(Rrs1_2+Rr_s_d+Ra_r1_2(1, i)/2)+(q2-p1)*(Ra_r1(1, i)+Ra_r1_2(1, i)/2)+q2*(2*Rrs1+2*Rr_s+Ra_r1(1, i))-(H1-D1/2)==0,...
        (H1-D1/2)-q1*(2*Rrs1+2*Rr_s+Ra_r1(1, i))-(q1+p1)*(Ra_r1(1, i)+Ra_s1/2)-(q1+p1+p3)*(Ra_s1/2+Ra_b1)-Hleaf==0];

Q=solve(eqns);
        obj.q1(1, i)=double(Q.q1); obj.q2(1, i)=double(Q.q2); obj.q3(1, i)=double(Q.q3);
        obj.p1(1, i)=double(Q.p1); obj.p2(1, i)=double(Q.p2); obj.p3(1, i)=double(Q.p3);

i=i+1
end

end


