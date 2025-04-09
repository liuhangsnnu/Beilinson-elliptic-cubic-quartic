/******************************************************
***********non-torsiont with 3 element over cubic field
*******************************************************/
{
\\ y^2=x^3+l^2*((p^2+p+1)*x+4*l*p^2*(p+1)^2)^2
\\ p+1 is a root of u^3+a*u^2+(a-1)*u+1
getcurvenontorsion(a,l) =
my(pol,rpol,q,u,sf2,sf3,K,L,K1,L1,M,M1,EK,EL);
pol=p^3+(a+3)*p^2+(a+2)*p+1;
K = nfinit(pol);
if(l!=4,
    rpol=rnfequation(K,x^2+4*l-l^2);
    L = nfinit(rpol);
    M = nfinit(x^2+4*l-l^2);
    sf2=nfsubfields(L,2);
    sf3=nfsubfields(L,3);
    K1 = nfinit(sf3[1][1]);
    M1 = nfinit(sf2[1][1]);
    \\ express p by x in L
    q=subst(nfisisom(K,K1)[1],x,sf3[1][2]);
    \\ express sqrt(l^2-4*l) by x in L
    u=subst(nfisisom(M,M1)[1],x,sf2[1][2]);
    \\ y^2=x^3+l^2*((p^2+p+1)*x+4*l*p^2*(p+1)^2)^2
    EL = ellinit([0,l^2*(q^2+q+1)^2,0,8*l^3*(q^2+q+1)*q^2*(q+1)^2,16*l^4*q^4*(q+1)^4], L);
);
EK = ellinit([0,l^2*(p^2+p+1)^2,0,8*l^3*(p^2+p+1)*p^2*(p+1)^2,16*l^4*p^4*(p+1)^4], K);
print("conductor norm of EK is ", factor(idealnorm(K, ellglobalred(EK)[1])));
return([EK,EL,q,u,K,L,idealnorm(K, ellglobalred(EK)[1]),K.disc]);
}

{
\\ compare regulator and L-value for the family with elements not supporting on torsion points
compareregLnontorsion(a,l) =
my(E,cond,EK,EL,K,L,q,u,d,OO,P1,P2,Q1,Q2,R1,R2,S1,S2,x1,x2,x3,tauK,l1,l2,l3,m,lstar,Clvec=[6,4,3,2],Cl=Clvec[l]);
E=getcurvenontorsion(a,l);
if(E==0, return);
EK=E[1];EL=E[2];q=E[3];u=E[4];K=E[5];L=E[6];cond=E[7];disc=E[8];
if(l==4,EL=EK;q=p;u=0);
\\ y^2=x^3+l^2*((p^2+p+1)*x+4*l*p^2*(p+1)^2)^2
\\ 3 torsion
P1 = [0,4*l^2*q^2*(q+1)^2];
\\ P2 = 2*P1
P2 = [0,-4*l^2*q^2*(q+1)^2];
\\ non 2 torsion
Q1 = [-q^2*4*l,4*q^3*l*u];
Q2 = [-q^2*4*l,-4*q^3*l*u];
R1 = [-(q+1)^2*4*l,4*(q+1)^3*l*u];
R2 = [-(q+1)^2*4*l,-4*(q+1)^3*l*u];
S1 = [-q^2*(q+1)^2*4*l,4*q^3*(q+1)^3*l*u];
S2 = [-q^2*(q+1)^2*4*l,-4*q^3*(q+1)^3*l*u];

tauK=calctau(EK);
tauL=calctau(EL);
print("tauK = ", tauK);
print("tauL = ", tauL);
xP=xLtoK(calcx(EL,P1),tauK,tauL);
x1=xLtoK(calcx(EL,elladd(EL,P1,Q1)),tauK,tauL);
x2=xLtoK(calcx(EL,elladd(EL,P1,S1)),tauK,tauL);
x3=xLtoK(calcx(EL,elladd(EL,P1,R1)),tauK,tauL);

x1m=xLtoK(calcx(EL,elladd(EL,P1,Q2)),tauK,tauL);
x2m=xLtoK(calcx(EL,elladd(EL,P1,S2)),tauK,tauL);
x3m=xLtoK(calcx(EL,elladd(EL,P1,R2)),tauK,tauL);

print("x1=",x1);
print("x2=",x2);
print("x3=",x3);

l1=6*Cl*(elldlog(x1,tauK)-I*elljlog(x1,tauK) + elldlog(x1m,tauK)-I*elljlog(x1m,tauK)-2*elldlog(xP,tauK)+2*I*elljlog(xP,tauK));
l2=6*Cl*(elldlog(x2,tauK)-I*elljlog(x2,tauK) + elldlog(x2m,tauK)-I*elljlog(x2m,tauK)-2*elldlog(xP,tauK)+2*I*elljlog(xP,tauK));
l3=6*Cl*(elldlog(x3,tauK)-I*elljlog(x3,tauK) + elldlog(x3m,tauK)-I*elljlog(x3m,tauK)-2*elldlog(xP,tauK)+2*I*elljlog(xP,tauK));


if(#tauK==3,
  \\ totally real
  m=[l1[1],l2[1],l3[1];l1[2],l2[2],l3[2];l1[3],l2[3],l3[3]],
  \\ not totally real
  ytau = imag(tauK[2]);
  m=[l1[1],l2[1],l3[1];\
  imag((l1[2]-conj(l1[2]))/ytau),imag((l2[2]-conj(l2[2]))/ytau),imag((l3[2]-conj(l3[2]))/ytau);\
  imag((tauK[2]*l1[2]-conj(tauK[2]*l1[2]))/ytau),imag((tauK[2]*l2[2]-conj(tauK[2]*l2[2]))/ytau),imag((tauK[2]*l3[2]-conj(tauK[2]*l3[2]))/ytau)]
);

my(lvalue=lfun(EK,0,3));
print("m=",matdet(m));
print("lvalue=",lvalue);
lstar=lvalue/3!;
dep=lindep([abs(2^numofrealtauhalf(tauK)*real(matdet(m))/(2*Pi)^3),lstar]);
return(totex(a,disc,cond,lstar,2*Cl*-dep[1]/dep[2]));

}

{
xLtoK(xL,tauK,tauL) =
my(xK=vector(#tauK));
for(i=1,#tauK,
  for(j=1,#tauL,
    if(abs(tauL[j]-tauK[i])<10^-10, xK[i]=xL[j]);
  )
);
return(xK);
}

/******************************************************
***********Z10-torsiont with 4 element over quatic field
*******************************************************/
{
getcurveZ10(a) =
my(pol,K,E,b,c,p);
\\ pol=u^4+a*u^3+14*u^2-a*u+1;
pol=u^4+a*u^3-a*u+1;
if(#factor(pol)[,1]>1,return(0));
print(pol);
K = nfinit(pol);
\\ u,(u+1)/(u-1) are units
\\ 10 torsion at [0,0]
p=(1-u)/2;
b=(2*p^5-3*p^4+p^3)/(p^2-3*p+1)^2;
c=(-2*p^3+3*p^2-p)/(p^2-3*p+1);
E = ellinit([(2*p^3 - 2*p^2 - 2*p + 1)*4,-p^3*(p-1)*(2*p-1)*4^2,-p^3*(p-1)*(2*p-1)*(p^2 - 3*p + 1)*4^3,0,0],K);
\\ 10 torsion point
P = [0,0];
print("conductor norm of E is ", factor(idealnorm(K, ellglobalred(E)[1])));
return([E,P,idealnorm(K, ellglobalred(E)[1]),K.disc]);
}

{
\\ compare regulator and L function for the family with Z/10 torsion
compareregLZ10(a) =
my(EP,cond,E,P1,P2,P3,P4,x1,x2,x3,x4,tau,ytau1,ytau2,l1,l2,l3,l4,m,lstar);
EP = getcurveZ10(a);
if(EP==0, return(""));
E = EP[1];
cond = EP[3];
disc = EP[4];
if(E==0, return);
\\ 10 torsion
P1 = EP[2];
P2 = elladd(E,P1,P1);
P3 = elladd(E,P1,P2);
P4 = elladd(E,P2,P2);

tau=calctau(E);
x1=calcx(E,P1);
x2=calcx(E,P2);
x3=calcx(E,P3);
x4=calcx(E,P4);

l1=2*10^3*(elldlog(x1,tau)-I*elljlog(x1,tau));
l2=2*10^3*(elldlog(x2,tau)-I*elljlog(x2,tau));
l3=2*10^3*(elldlog(x3,tau)-I*elljlog(x3,tau));
l4=2*10^3*(elldlog(x4,tau)-I*elljlog(x4,tau));

if(#tau==4,
  \\ totally real
  m=[l1[1],l2[1],l3[1],l4[1];l1[2],l2[2],l3[2],l4[2];l1[3],l2[3],l3[3],l4[3];l1[4],l2[4],l3[4],l4[4]],
);
if(#tau==3,  
  \\ two real embeddings and one complex embedding
  ytau = imag(tau[3]);
  m=[l1[1],l2[1],l3[1],l4[1];\
  l1[2],l2[2],l3[2],l4[2];\
  2*imag(l1[3])/ytau,2*imag(l2[3])/ytau,2*imag(l3[3])/ytau,2*imag(l4[3])/ytau;\
  2*real(l1[3]),2*real(l2[3]),2*real(l3[3]),2*real(l4[3])]
);
if(#tau==2,  
  \\ two complex embeddings
  ytau1 = imag(tau[1]);
  ytau2 = imag(tau[2]);
  m=[2*imag(l1[1])/ytau1,2*imag(l2[1])/ytau1,2*imag(l3[1])/ytau1,2*imag(l4[1])/ytau1;\
  2*real(l1[1]),2*real(l2[1]),2*real(l3[1]),2*real(l4[1]);\
  2*imag(l1[2])/ytau2,2*imag(l2[2])/ytau2,2*imag(l3[2])/ytau2,2*imag(l4[2])/ytau2;\
  2*real(l1[2]),2*real(l2[2]),2*real(l3[2]),2*real(l4[2])]
);

print("tau=",tau);
print("x1=",x1);
print("x2=",x2);
print("x3=",x3);
print("x4=",x4);
print("l1=",l1);
print("l2=",l2);
print("l3=",l3);
print("l4=",l4);
print(matdet(m));


my(lvalue=lfun(E,0,4));
print(lvalue);
lstar=lvalue/4!;
dep=lindep([abs(2^numofrealtauhalf(tau)*real(matdet(m))/(2*Pi)^4),lstar]);
return(totex(a,disc,cond,lstar,-dep[1]/dep[2]));

}

/******************************************************
***********Z8-torsiont with 3 element over cubic field
*******************************************************/
{
compareregLZ8(a) =
\\ compare regulator and L function for the family with Z/8 torsion
my(E,cond,P1,P2,P3,x1,x2,x3,tau,l1,l2,l3,m,lstar);
[E,cond,disc]=getcurveZ8(a);
if(E==0, return);
\\ 8 torsion
P1 = [0,0];
P2 = elladd(E,P1,P1);
P3 = elladd(E,P2,P1);
print("P1=", P1);
print("P2=", P2);
print("P3=", P3);
tau=calctau(E);
x1=calcx(E,P1);
x2=calcx(E,P2);
x3=calcx(E,P3);
l1=2*8^3*(elldlog(x1,tau)-I*elljlog(x1,tau));
l2=2*8^3*(elldlog(x2,tau)-I*elljlog(x2,tau));
l3=2*8^3*(elldlog(x3,tau)-I*elljlog(x3,tau));

print("tau=",tau);
print("x1=",x1);
print("x2=",x2);
print("x3=",x3);
print("l1=",l1);
print("l2=",l2);
print("l3=",l3);

if(#tau==3,
  \\ totally real
  m=[l1[1],l2[1],l3[1];l1[2],l2[2],l3[2];l1[3],l2[3],l3[3]],
  \\ not totally real
  ytau = imag(tau[2]);
  m=[l1[1],l2[1],l3[1];\
  2*imag(l1[2])/ytau,2*imag(l2[2])/ytau,2*imag(l3[2])/ytau;\
  2*real(l1[2]),2*real(l2[2]),2*real(l3[2])]
);

my(lvalue=lfun(E,0,3));
print("m=",matdet(m));
print("lvalue=",lvalue);
lstar=lvalue/3!;
dep=lindep([abs(2^numofrealtauhalf(tau)*real(matdet(m))/(2*Pi)^3),lstar]);
return(totex(a,disc,cond,lstar,-dep[1]/dep[2]));
}


{
getcurveZ8(a) =
my(pol,K,E,n,d);
pol=p^3+a*p^2+(-a-1)*p+1;
if(#factor(pol)[,1]>1,return(0));
K = nfinit(pol);
t=1/Mod(p+1,pol);
b=2*t^2-3*t+1;
c=(2*t^2-3*t+1)*(p+1);
E=ellinit([1-c,-b,-b,0,0],K);
E=ellchangecurve(E,[t,0,0,0]);
print("conductor norm of E is ", factor(idealnorm(K, ellglobalred(E)[1])));
return([E,idealnorm(K, ellglobalred(E)[1]),K.disc]);
}





/******************************************************
***********7-torsiont with 3 element over cubic field
*******************************************************/
{
compareregLZ7(a) =
\\ compare regulator and L function for the family with Z/7 torsion 
my(E,cond,P1,P2,P3,x1,x2,x3,tau,l1,l2,l3,m,lstar);
[E,cond,disc]=getcurveZ7(a);
if(E==0, return);
\\ 7 torsion
P1 = [0,0];
P2 = elladd(E,P1,P1);
P3 = elladd(E,P2,P1);
print("P1=", P1);
print("P2=", P2);
print("P3=", P3);
tau=calctau(E);
x1=calcx(E,P1);
x2=calcx(E,P2);
x3=calcx(E,P3);
l1=7^3*(elldlog(x1,tau)-I*elljlog(x1,tau));
l2=7^3*(elldlog(x2,tau)-I*elljlog(x2,tau));
l3=7^3*(elldlog(x3,tau)-I*elljlog(x3,tau));

print("tau=",tau);
print("x1=",x1);
print("x2=",x2);
print("x3=",x3);
print("l1=",l1);
print("l2=",l2);
print("l3=",l3);
if(#tau==3,
  \\ totally real
  m=[l1[1],l2[1],l3[1];l1[2],l2[2],l3[2];l1[3],l2[3],l3[3]],
  \\ not totally real
  ytau = imag(tau[2]);

  m=[l1[1],l2[1],l3[1];\
  2*imag(l1[2])/ytau,2*imag(l2[2])/ytau,2*imag(l3[2])/ytau;\
  2*real(l1[2]),2*real(l2[2]),2*real(l3[2])]
);

my(lvalue=lfun(E,0,3));
print("m=",matdet(m));
print("lvalue=",lvalue);
\\ for real embedding with real(tau) = 1/2, the generator of H^- is 2*tau-1, multiply 2
lstar=lvalue/3!;
dep=lindep([abs(2^numofrealtauhalf(tau)*real(matdet(m))/(2*Pi)^3),lstar]);
return(totex(a,disc,cond,lstar,-dep[1]/dep[2]));
}

{
getcurveZ7(a) =
my(pol,K,E,n,d);
pol=p^3+a*p^2+(-a-1)*p+1;
if(#factor(pol)[,1]>1,return(0));
K = nfinit(pol);
b=p^3-p^2;
c=p^2-p;
E=ellinit([1-c,-b,-b,0,0],K);
print("conductor norm of E is ", factor(idealnorm(K, ellglobalred(E)[1])));
return([E, idealnorm(K, ellglobalred(E)[1]), K.disc]);
}

{
numofrealtauhalf(tau) =
my(num=0);
for(i=1, #tau, if(abs(frac(real(tau[i]))-1/2)<10^-10, num=num+1));
return(num);
}

/******************************************************
***********calculate elliptic dilogrithm
*******************************************************/
{
elldlog(x,tau)  =
if(#x!=#tau,print("length of x and tau are not the same"));
my(result=vector(#tau));
for(i=1,#tau,
  my(q=exp(2*Pi*I*tau[i]));
  my(z=exp(2*Pi*I*x[i]));
  for(n=-1000,1000,result[i] += polylog(2,z*q^n,2));
);
return(result);
}

{
elljlog(x,tau)  =
if(#x!=#tau,print("length of x and tau are not the same"));
my(result=vector(#tau));
for(i=1,#tau,
  my(q=exp(2*Pi*I*tau[i]));
  my(z=exp(2*Pi*I*x[i]));
  for(n=0,1000,result[i] += j(z*q^n));
  for(n=1,1000,result[i] -= j(z^(-1)*q^n));
  result[i] += 1/3*log(abs(q))^2*subst(bernpol(3,y),y,log(abs(z))/log(abs(q)));
);
return(result);
}


{
j(z) =
return(log(abs(z))*log(abs(1-z)));
}

{
\\ E defined over real number
calctau(E) =
my(o=E.omega);
my(tau);
if(type(o[1])=="t_REAL",tau=vector(1);tau[1]=o[2]/o[1];
  if(imag(tau[1])<0,tau[1]=-tau[1]);
);
if(type(o[1])=="t_VEC",
  tau=vector(#o);
  for(i=1,#o,
    tau[i]=o[i][2]/o[i][1];
    \\ if the o[i][1] is imaginary, this happens for getcurvenontorsion
    if(abs(imag(o[i][1]))>1E-10,tau[i]=o[i][1]/o[i][2]);
    if(imag(tau[i])<0,tau[i]=-tau[i]);
  );
);
return(tau);
}

{
calcx(E, pt) =
my(o=E.omega);
my(x);
if(type(o[1])=="t_REAL",x=vector(1);x[1]=ellpointtoz(E,pt)/o[1]);
if(type(o[1])=="t_VEC",
  x=vector(#o);
  for(i=1,#o,
    x[i]=ellpointtoz(E,pt)[i]/o[i][1];
    \\ if the o[i][1] is imaginary, this happens for getcurvenontorsion
    if(abs(imag(o[i][1]))>1E-10,x[i]=ellpointtoz(E,pt)[i]/o[i][2];);
  );
);
return(x);
}

{
totex(a, disc, cond, Lstar, Qa) =
my(s="");
s=Str("$",a,"$ & ", totexint(disc), " & ", totexint(cond), " & $", Lstar, "$ & ", totexint(Qa), "\\\\");
}

{
totexint(c) =
my(fac, s="$");
fac=factor(c);
ind=1;
if(fac[1,1]==-1,s=Str(s,"-");ind=2);
for(i=ind,#fac[,1],
s=Str(s, fac[i,1]);
if(fac[i,2]!=1,s=Str(s, "^{",fac[i,2],"}"));
if(i!=#fac[,1],s=Str(s, " \\cdot "))
);
s=Str(s,"$");
return(s);
}

{
gettable7(abegin, aend) =
my(s="");
for(i=abegin,aend,
s=Str(s,"\n",compareregLZ7(i))
);
print(s);
}

{
gettable8(abegin, aend) =
my(s="");
for(i=abegin,aend,
s=Str(s,"\n",compareregLZ8(i))
);
print(s);
}

{
gettable10(abegin, aend) =
my(s="");
for(i=abegin,aend,
s=Str(s,"\n",compareregLZ10(i))
);
print(s);
}

{
gettablenontorsion(abegin, aend, l) =
my(s="");
for(i=abegin,aend,
s=Str(s,"\n",compareregLnontorsion(i,l))
);
print(s);
}

