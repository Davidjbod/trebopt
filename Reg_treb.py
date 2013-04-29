#Floating Arm Trebuchet Trim   David Bodkin  v 0.1
#Equations from: http://students.eou.edu/~duyckj/treb-final/node4.html

#Current assumptions
#Frictionless, point masses, wheel diameter represented by single point.

import sys 

from numpy import *
import matplotlib.pyplot as plt
import time as time

#import cProfile
#import pstats




class opts:
    def __init__(self, opt1, opt2):
        self.interactive = opt1
        self.drag = opt2



def dth2(m1,m2,mb,l1,l2,l3,l4,l5,Lb,th,phi,psi,psis,thd,phid,psid,flag):
    sth=sin(th)
    cth=cos(th)
    sphi=sin(phi)
    cphi=cos(phi)
    spsi=sin(psi)
    cpsi=cos(psi)
    
    spsis=sin(psis)
    cpsis=cos(psis)
    
    spsissq=spsis**2

    thdsq=thd*thd
    phidsq=phid*phid

    sthsq=sth*sth
    sthcu=sthsq*sth
    cthsq=cth*cth
    cthcu=cthsq*cth
    cphisq=cphi*cphi

    l1sq=l1**2
    l1cu=l1**3
    l1qt=l1**4
    l1qn=l1**5
    
    l2sq=l2**2
    l2cu=l2**3
    l2qt=l2**4
    l2qn=l2**5

    l3sq=l3**2
    l3cu=l3**3
    l3qt=l3**4
    l3qn=l3**5
    
    if flag == 0:  #Projectile on  ground
        
        rad= sqrt(1-(l2**2*(cth + spsis)*(cth + spsis))/l3**2)
        r2=rad*rad
        r3=r2*rad
        r4=r2*r2

        # for th in the constrained portion of the throw:
        c1=-6*l2cu*l3sq*m2*r3*thdsq*cthcu + 6*l1*l3cu*l4*m1*phidsq*r4*sphi
        c2= 12*l1*l3cu*l4*m1*phid*r4*thd*sphi + 6*l1*l3cu*l4*m1*r4*thdsq*sphi
        c3=-6*l1sq*l3cu*m1*r4*thdsq*cphi*sphi
        c4=-6*l2cu*l3sq*m2*r3*thdsq*cthsq*spsis + 6*g*l1*l3cu*m1*r4*sth
        c5= 3*g*l1*l3cu*mb*r4*sth-3*g*l2*l3cu*mb*r4*sth
        c6= 6*l2sq*l3cu*m2*r2*thdsq*cth*sth
        c7=-12*l2sq*l3cu*m2*r4*thdsq*cth*sth
        c8= 12*l2cu*l3sq*m2*r3*thdsq*cth*sthsq
        c9= 6*l2qn*m2*rad*thdsq*cthcu*sthsq
        c10= 6*l2cu*l3sq*m2*r3*thdsq*spsis*sthsq
        c11= 12*l2qn*m2*rad*thdsq*cthsq*spsis*sthsq
        c12= 6*l2qn*m2*rad*thdsq*cth*spsissq*sthsq
        c13=-6*l2qt*l3*m2*thdsq*cth*sthcu-6*l2qt*l3*m2*thdsq*spsis*sthcu
        c14=-6*g*l1*l3cu*m1*r4*cphi*sin(phi + th)
        cnum=c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14

        d1=-6*l1sq*l3cu*m1*r4-6*l2sq*l3cu*m2*r4-2*l1sq*l3cu*mb*r4
        d2= 2*l1*l2*l3cu*mb*r4-2*l2sq*l3cu*mb*r4 + 6*l1sq*l3cu*m1*r4*cphisq 
        d3= 12*l2cu*l3sq*m2*r3*cthsq*sth
        d4= 12*l2cu*l3sq*m2*r3*cth*spsis*sth-6*l2sq*l3cu*m2*r2*sthsq
        d5= 12*l2sq*l3cu*m2*r4*sthsq
        ddenom=d1+d2+d3+d4+d5

        thdd=cnum/ddenom
    else:	#In air
        a1=6*l1*l4*m1*phid*phid*sphi + 12*l1*l4*m1*phid*thd*sphi
        a2=6*l1*l4*m1*thd*thd*sphi-6*l1sq*m1*thd*thd*cphi*sphi
        a3=-6*l2*l3*m2*psid*psid*spsi + 12*l2*l3*m2*psid*thd*spsi
        a4=-6*l2*l3*m2*thd*thd*spsi + 6*l2sq*m2*thd*thd*cpsi*spsi
        a5=-6*g*l2*m2*cpsi*(spsi*cth-sth*cpsi) + 6*g*l1*m1*sth-6*g*l2*m2*sth
        a6=3*g*l1*mb*sth-3*g*l2*mb*sth-6*g*l1*m1*cphi*(sphi*cth+sth*cphi)
        thtop=a1+a2+a3+a4+a5+a6 

        b1=-6*l1sq*m1-6*l2sq*m2-2*l1sq*mb + 2*l1*l2*mb-2*l2sq*mb
        b2=6*l1sq*m1*cphi*cphi + 6*l2sq*m2*cpsi*cpsi
        thbot=b1+b2
        thdd=thtop/thbot

    return thdd


def dphi2(m1,m2,mb,l1,l2,l3,l4,l5,Lb,th,phi,psi,psis,thd,phid,psid,flag)  :
    

    sth=sin(th)
    cth=cos(th)
    sphi=sin(phi)
    cphi=cos(phi)
    spsi=sin(psi)
    cpsi=cos(psi)

    spsis=sin(psis)
    cpsis=cos(psis)
    spsissq=spsis**2

    thdsq=thd*thd
    phidsq=phid*phid
    sthsq=sth*sth
    sthcu=sthsq*sth
    cthsq=cth*cth
    cthcu=cthsq*cth
    cphisq=cphi*cphi

    l1sq=l1**2
    l1cu=l1**3
    l1qt=l1**4
    l1qn=l1**5
    
    l2sq=l2**2
    l2cu=l2**3
    l2qt=l2**4
    l2qn=l2**5

    l3sq=l3**2
    l3cu=l3**3
    l3qt=l3**4
    l3qn=l3**5

    l4sq=l4**2
    l4cu=l4**3
    l4qt=l4**4
    l4qn=l4**5

    m1sq=m1**2;

    a=zeros([35,1])
    b=zeros([5,1])
    if flag==0:
       

        rad= sqrt(1-(l2**2*(cth + spsis)*(cth + spsis))/l3**2)
        r2=rad*rad
        r3=r2*rad
        r4=r2*r2

        a[1]=-(l1*l4cu*m1sq*phidsq*sphi) + cphi*l1sq*l4sq*m1sq*phidsq*sphi +cphi*g*l1sq*l4*m1sq*sth
        a[2]=-g*l1*l4sq*m1sq*sth + (cphi*g*l1sq*l4*m1*mb*sth)/2.-(cphi*g*l1*l2*l4*m1*mb*sth)/2.
        a[3]=-(g*l1*l4sq*m1*mb*sth)/2.+ (g*l2*l4sq*m1*mb*sth)/2.-2*l1*l4cu*m1sq*phid*sphi*thd
        a[4]= 2*cphi*l1sq*l4sq*m1sq*phid*sphi*thd-(cphi*cthcu*l1*l2cu*l4*m1*m2*thdsq)/(l3*rad)
        a[5]= (cthcu*l2cu*l4sq*m1*m2*thdsq)/(l3*rad)-l1cu*l4*m1sq*sphi*thdsq
        a[6]=-l1*l4cu*m1sq*sphi*thdsq + 2*cphi*l1sq*l4sq*m1sq*sphi*thdsq-l1*l2sq*l4*m1*m2*sphi*thdsq
        a[7]=-(l1cu*l4*m1*mb*sphi*thdsq)/3. + (l1sq*l2*l4*m1*mb*sphi*thdsq)/3.
        a[8]=-(l1*l2sq*l4*m1*mb*sphi*thdsq)/3.-(cphi*cthsq*l1*l2cu*l4*m1*m2*spsis*thdsq)/(l3*rad)
        a[9]= (cthsq*l2cu*l4sq*m1*m2*spsis*thdsq)/(l3*rad)-2*cphi*cth*l1*l2sq*l4*m1*m2*sth*thdsq
        a[10]= 2*cth*l2sq*l4sq*m1*m2*sth*thdsq +(cphi*cth*l1*l2sq*l4*m1*m2*sth*thdsq)/r2
        a[11]=-(cth*l2sq*l4sq*m1*m2*sth*thdsq)/r2
        a[12]= (2*cthsq*l1*l2cu*l4*m1*m2*sphi*sth*thdsq)/(l3*rad)
        a[13]= (2*cth*l1*l2cu*l4*m1*m2*sphi*spsis*sth*thdsq)/(l3*rad)
        a[14]=-(cphi*cth*l1*l2qt*l4*m1*m2*sthcu*thdsq)/(l3sq*r4)
        a[15]= (cth*l2qt*l4sq*m1*m2*sthcu*thdsq)/(l3sq*r4)
        a[16]=-(cphi*l1*l2qt*l4*m1*m2*spsis*sthcu*thdsq)/(l3sq*r4)
        a[17]= (l2qt*l4sq*m1*m2*spsis*sthcu*thdsq)/(l3sq*r4)
        a[18]= (cphi*cthcu*l1*l2qn*l4*m1*m2*sthsq*thdsq)/(l3cu*r3)
        a[19]=-(cthcu*l2qn*l4sq*m1*m2*sthsq*thdsq)/(l3cu*r3)
        a[20]= (2*cphi*cth*l1*l2cu*l4*m1*m2*sthsq*thdsq)/(l3*rad)
        a[21]=-(2*cth*l2cu*l4sq*m1*m2*sthsq*thdsq)/(l3*rad) +2*l1*l2sq*l4*m1*m2*sphi*sthsq*thdsq
        a[22]=-(l1*l2sq*l4*m1*m2*sphi*sthsq*thdsq)/r2
        a[23]= (2*cphi*cthsq*l1*l2qn*l4*m1*m2*spsis*sthsq*thdsq)/(l3cu*r3) 
        a[24]=-(2*cthsq*l2qn*l4sq*m1*m2*spsis*sthsq*thdsq)/(l3cu*r3)
        a[25]= (cphi*l1*l2cu*l4*m1*m2*spsis*sthsq*thdsq)/(l3*rad)
        a[26]=-(l2cu*l4sq*m1*m2*spsis*sthsq*thdsq)/(l3*rad)
        a[27]= (cphi*cth*l1*l2qn*l4*m1*m2*spsissq*sthsq*thdsq)/(l3cu*r3)
        a[28]=-(cth*l2qn*l4sq*m1*m2*spsissq*sthsq*thdsq)/(l3cu*r3)-g*l1sq*l4*m1sq*sin(phi + th)
        a[29]= cphi*g*l1*l4sq*m1sq*sin(phi + th)-g*l2sq*l4*m1*m2*sin(phi + th)
        a[30]=-(g*l1sq*l4*m1*mb*sin(phi + th))/3. + (g*l1*l2*l4*m1*mb*sin(phi +th))/3.
        a[31]=-(g*l2sq*l4*m1*mb*sin(phi + th))/3.
        a[32]= (2*cthsq*g*l2cu*l4*m1*m2*sth*sin(phi + th))/(l3*rad)
        a[33]= (2*cth*g*l2cu*l4*m1*m2*spsis*sth*sin(phi + th))/(l3*rad)
        a[34]= 2*g*l2sq*l4*m1*m2*sthsq*sin(phi + th)-(g*l2sq*l4*m1*m2*sthsq*sin(phi + th))/r2

        b[1]=-(l1sq*l4sq*m1sq) + cphisq*l1sq*l4sq*m1sq-l2sq*l4sq*m1*m2-(l1sq*l4sq*m1*mb)/3.
        b[2]= (l1*l2*l4sq*m1*mb)/3.-(l2sq*l4sq*m1*mb)/3. +(2*cthsq*l2cu*l4sq*m1*m2*sth)/(l3*rad)
        b[3]= (2*cth*l2cu*l4sq*m1*m2*spsis*sth)/(l3*rad) + 2*l2sq*l4sq*m1*m2*sthsq
        b[4]=-(l2sq*l4sq*m1*m2*sthsq)/r2
        phiddtopc=0.


        for i in range(1,35):
            phiddtopc=phiddtopc+a[i]
        next
        
        phiddbotc=b[1]+b[2]+b[3]+b[4]

        phidd=phiddtopc/phiddbotc
        

        
    else:
        c1=-6*l1*l4sq*m1*phid*phid*sphi-12*l1*l4sq*m1*phid*thd*sphi-6*l1cu*m1*thd*thd*sphi
        c2=-6*l1*l4sq*m1*thd*thd*sphi-6*l1*l2sq*m2*thd*thd*sphi-2*l1cu*mb*thd*thd*sphi
        c3=2*l1sq*l2*mb*thd*thd*sphi-2*l1*l2sq*mb*thd*thd*sphi +6*l1sq*l4*m1*phid*phid*cphi*sphi
        c4=12*l1sq*l4*m1*phid*thd*cphi*sphi + 12*l1sq*l4*m1*thd*thd*cphi*sphi
        c5=6*l1*l2sq*m2*thd*thd*cpsi*cpsi*sphi + 6*l2*l3*l4*m2*psid*psid*spsi
        c6=-12*l2*l3*l4*m2*psid*thd*spsi + 6*l2*l3*l4*m2*thd*thd*spsi
        c7=-6*l1*l2*l3*m2*psid*psid*cphi*spsi + 12*l1*l2*l3*m2*psid*thd*cphi*spsi
        c8=-6*l1*l2*l3*m2*thd*thd*cphi*spsi-6*l2sq*l4*m2*thd*thd*cpsi*spsi
        c9=6*l1*l2sq*m2*thd*thd*cphi*cpsi*spsi + 6*g*l2*l4*m2*cpsi*(spsi*cth-sth*cpsi)
        c10=-6*g*l1*l2*m2*cphi*cpsi*(spsi*cth-sth*cpsi)-6*g*l1*l4*m1*sth +6*g*l2*l4*m2*sth
        c11=-3*g*l1*l4*mb*sth + 3*g*l2*l4*mb*sth + 6*g*l1sq*m1*cphi*sth
        c12=-6*g*l1*l2*m2*cphi*sth + 3*g*l1sq*mb*cphi*sth-3*g*l1*l2*mb*cphi*sth
        c13=-6*g*l1sq*m1*(sphi*cth+sth*cphi)-6*g*l2sq*m2*(sphi*cth+sth*cphi)-2*g*l1sq*mb*(sphi*cth+sth*cphi)
        c14=2*g*l1*l2*mb*(sphi*cth+sth*cphi)-2*g*l2sq*mb*(sphi*cth+sth*cphi) +6*g*l1*l4*m1*cphi*(sphi*cth+sth*cphi)
        c15=6*g*l2sq*m2*cpsi*cpsi*(sphi*cth+sth*cphi)

        d1=-6*l1sq*l4*m1-6*l2sq*l4*m2-2*l1sq*l4*mb + 2*l1*l2*l4*mb-2*l2sq*l4*mb
        d2=6*l1sq*l4*m1*cphi*cphi + 6*l2sq*l4*m2*cpsi*cpsi

        phitop=c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14+c15
        phibot=d1+d2
        phidd=phitop/phibot

    return phidd



def dpsi2(m1,m2,mb,l1,l2,l3,l4,l5,Lb,th,phi,psi,psis,thd,phid,psid,flag)  :
    

    sth=sin(th)
    cth=cos(th)
    sphi=sin(phi)
    cphi=cos(phi)
    spsi=sin(psi)
    cpsi=cos(psi)

    thdsq=thd*thd
    phidsq=phid*phid
    sthsq=sth*sth
    sthcu=sthsq*sth
    cthsq=cth*cth
    cthcu=cthsq*cth
    cphisq=cphi*cphi

    l1sq=l1**2
    l1cu=l1**3
    l1qt=l1**4
    l1qn=l1**5
    
    l2sq=l2**2
    l2cu=l2**3
    l2qt=l2**4
    l2qn=l2**5

    l3sq=l3**2
    l3cu=l3**3
    l3qt=l3**4
    l3qn=l3**5

    l4sq=l4**2
    l4cu=l4**3
    l4qt=l4**4
    l4qn=l4**5

    m1sq=m1**2
    m2sq=m2**2;

    if flag==0:
        psidd=0
    else:
        e1=l1*l3sq*l4cu*m1sq*m2*phid*phid*sphi + 2*l1*l3sq*l4cu*m1sq*m2*phid*thd*sphi
        e2=l1*l3sq*l4cu*m1sq*m2*thd*thd*sphi-l1sq*l3sq*l4sq*m1sq*m2*thd*thd*cphi*sphi
        e3=-l1*l2*l3*l4cu*m1sq*m2*phid*phid*cpsi*sphi-2*l1*l2*l3*l4cu*m1sq*m2*phid*thd*cpsi*sphi
        e4=-l1*l2*l3*l4cu*m1sq*m2*thd*thd*cpsi*sphi
        e5=l1sq*l2*l3*l4sq*m1sq*m2*thd*thd*cphi*cpsi*sphi-l2*l3cu*l4sq*m1*m2sq*psid*psid*spsi
        e6=2*l2*l3cu*l4sq*m1*m2sq*psid*thd*spsi-l1sq*l2*l3*l4sq*m1sq*m2*thd*thd*spsi
        e7=-l2cu*l3*l4sq*m1*m2sq*thd*thd*spsi-l2*l3cu*l4sq*m1*m2sq*thd*thd*spsi
        e8=-(l1sq*l2*l3*l4sq*m1*m2*mb*thd*thd*spsi)/3. +(l1*l2sq*l3*l4sq*m1*m2*mb*thd*thd*spsi)/3.
        e9=-(l2cu*l3*l4sq*m1*m2*mb*thd*thd*spsi)/3. +l1sq*l2*l3*l4sq*m1sq*m2*thd*thd*cphi*cphi*spsi
        e10=l2sq*l3sq*l4sq*m1*m2sq*psid*psid*cpsi*spsi
        e11=-2*l2sq*l3sq*l4sq*m1*m2sq*psid*thd*cpsi*spsi
        e12=2*l2sq*l3sq*l4sq*m1*m2sq*thd*thd*cpsi*spsi +g*l1sq*l3*l4sq*m1sq*m2*(spsi*cth-sth*cpsi)
        e13=g*l2sq*l3*l4sq*m1*m2sq*(spsi*cth-sth*cpsi) +(g*l1sq*l3*l4sq*m1*m2*mb*(spsi*cth-sth*cpsi))/3.
        e14=-(g*l1*l2*l3*l4sq*m1*m2*mb*(spsi*cth-sth*cpsi))/3. +(g*l2sq*l3*l4sq*m1*m2*mb*(spsi*cth-sth*cpsi))/3.
        e15=-g*l1sq*l3*l4sq*m1sq*m2*cphi*cphi*(spsi*cth-sth*cpsi)-g*l2*l3sq*l4sq*m1*m2sq*cpsi*(spsi*cth-sth*cpsi)
        e16=g*l1*l3sq*l4sq*m1sq*m2*sth-g*l2*l3sq*l4sq*m1*m2sq*sth
        e17=(g*l1*l3sq*l4sq*m1*m2*mb*sth)/2.-(g*l2*l3sq*l4sq*m1*m2*mb*sth)/2.
        e18=-g*l1*l2*l3*l4sq*m1sq*m2*cpsi*sth + g*l2sq*l3*l4sq*m1*m2sq*cpsi*sth
        e19=-(g*l1*l2*l3*l4sq*m1*m2*mb*cpsi*sth)/2. +(g*l2sq*l3*l4sq*m1*m2*mb*cpsi*sth)/2.
        e20=-g*l1*l3sq*l4sq*m1sq*m2*cphi*(sphi*cth+sth*cphi) +g*l1*l2*l3*l4sq*m1sq*m2*cphi*cpsi*(sphi*cth+sth*cphi)

        psitop=e1+e2+e3+e4+e5+e6+e7+e8+e9+e10+e11+e12+e13+e14+e15+e16+e17+e18+e19+e20
        f1=-(l1sq*l3sq*l4sq*m1sq*m2)-l2sq*l3sq*l4sq*m1*m2sq-(l1sq*l3sq*l4sq*m1*m2*mb)/3.
        f2=(l1*l2*l3sq*l4sq*m1*m2*mb)/3.-(l2sq*l3sq*l4sq*m1*m2*mb)/3.
        f3=l1sq*l3sq*l4sq*m1sq*m2*cphi*cphi
        f4=l2sq*l3sq*l4sq*m1*m2sq*cpsi*cpsi
        
        psibot=f1+f2+f3+f4
        
        psidd=psitop/psibot

    return psidd




#Find max val and location of max val
def maxval(arr):
    val=-10**-12
    index=0
    for i in range(0,len(arr)-1):
        if arr[i]>val:
            val=arr[i]
            index=i

    return [val,index]


def accelx(Vx):
    out=-(0.5*rho*Vx**2*S*cd)/m
    return out

def accely(Vy):
    out=-(0.5*rho*Vy**2*S*cd)/m*Vy/abs(Vy)-g
    return out

#Get input
def getinputfromfile():
    

    #Trebuchet Setup
    if len(sys.argv)>=2:
        file=sys.argv[1]
    else:
        print '***Error! - File Name required in command line'
        sys.exit()

    try:
        infile = open(file,'r')
    except:
        print '***Error! - Invalid file name'
        sys.exit()
        
    indata=[]
    for line in infile:
        if line != '\n':
            indata.append(line.split()[0])
            
    
    g=32.17;      #Gravity - ft/s^2
    M=float(indata[6])/g;      #Counter weight mass - slugs
    m=float(indata[7])/g;        #Projectile mass - slugs

    l1=float(indata[2])	      #Distance from tip of throwing arm to axle - ft
    l2=float(indata[1])	      #Distance from counter weight cg to axle - ft
    l3=float(indata[3])            #Sling Length - ft
    l5=float(indata[4])              #Axle Height - ft
    l4=float(indata[5])              #CW arm length - ft

    mb=float(indata[8])/g         #Bar mass
    if float(indata[9])==-1.0:
        Lb=(l1+l2)/3.             #Bar CG
    else:
        Lb=float(indata[9])

    #Sim Setup
    hmax=float(indata[13])       #Max dt for Runge-Kutta Method
    hmin=float(indata[14])       #Min dt for RK
    Tol=float(indata[15])        #Tolerance for R error
    N=int(indata[16])  	         #Max number of iterations

    runopts.interactive = int(indata[18])
    runopts.drag = int(indata[19])
    runopts.verbose = int(indata[20])
    runopts.gss=0

    #Inital values

    if float(indata[10])==999:
        thetai=arcsin(l5/l2)+pi/2
    else:    
        thetai=float(indata[10])*pi/180.0	#Aangle between arm and horizontal

    dia=float(indata[11])

    return [M,m,mb,l1,l2,l3,l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]

#Main code
def main(inputs,runopts):
    global rho
    global g
    global m
    global S
    global cd
    
    start = time.time()
    

    [M,m,mb,l1,l2,l3,l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]=inputs
    
    thetai=arcsin(l5/l2)+pi/2
    

    w=0               #Internal variable for RK for theta
    w2=0              #Internal variable for RK for phi
    w3=0              #Internal variable for RK for psi
    
    #Initialize arrays and matricies

    x3=zeros([N,2])	#projectile position in sling
    v3=zeros([N,2])	#projective velocity in sling
    a3=zeros([N,2])	#projective acceleration
    
    x4=zeros([N,2])	#Counter weight position
    v4=zeros([N,2])	#Counter weight velocity
    
    
    xtip=zeros([N,2])
    x1=zeros([N,2])
    xbarcg=zeros([N,2])

    dthdt2=zeros([N,1])	#Derivative of Theta wrt time
    dphidt2=zeros([N,1])	#Derivative of Phi wrt time
    dpsidt2=zeros([N,1])	#Derivative of Psi wrt time
    
    dthdt=zeros([N,1])	#Derivative of Theta wrt time
    dphidt=zeros([N,1])	#Derivative of Phi wrt time
    dpsidt=zeros([N,1])	#Derivative of Psi wrt time

    theta=zeros([N,1])	#Angle between arm and horizontal (0 deg = horizonal, Neg. angles = CW above axle
    phi=zeros([N,1])	#Angle between arm and CW arm (0deg = aligned with arm overlapping)
    psi=zeros([N,1])    #Angle arm and sling (0 deg = sling aligned with arm, Pos. angles = normal movement *180/pi
    
    T=zeros([N,1])	    #Time - s
    Ten=zeros([N,1])	#Sling Tension - lbf
    Ten2=zeros([N,1])	#Sling Tension - lbf
    Fb=zeros([N,1])	    #Force on projectile in the vertical direction (Pos. = up) - lbf
    Rnd=zeros([N,1])	#Calculated zero drag range if released - ft

    PEcw=zeros([N,1])	#Potential Energy of Counter weight
    PEp =zeros([N,1])	#Potential Energy of Projectile
    PEbar=zeros([N,1])	#Potenial Energy of bar    
    KEcw=zeros([N,1])	#Kenetic Energy of Counter weight
    KEp=zeros([N,1])	#Kenetic Energy of Projectile
    KEbar=zeros([N,1])	#Kenetic Energy of bar  
      
    KEtot=zeros([N,1])	#Kenetic Energy total  
    PEtot=zeros([N,1])	#Potential Energy total 
    
    Etot=zeros([N,1])   #Total Energy

    Rstore=zeros([N,1])     #RK Error diff
    hstore=zeros([N,1])     #h storage
    rkiterstore=zeros([N,1])     #h storage

    Fp=zeros([N,2])         #Force on bar due to sling tension
    Farmw=zeros([N,2])      #Force on bar due to bar's weight
    Fcw=zeros([N,2])        #Force on bar due to CW
    Rpin=zeros([N,2])        #Force on bar due to CW
    x=0                     #location for max moment in moment eqs
    #Fyp4=zeros([N,1])       #Bar force in yp due to Mb*g
    #Fyp5=zeros([N,1])       #Bar force in yp due to T

    Mombar1=zeros([N,1])
    Mombar2=zeros([N,1])
    Mombar3=zeros([N,1])

    angout=zeros([N,1])
    Rprime=zeros([N,1])

    out=zeros([1000,2])
    rvel=zeros([1000,2])
    Pos=zeros([1000,2])
    Rwd=zeros([N,1])
    Tp=zeros([1000,1])
    
  
    


    g=32.17;      #Gravity - ft/s^2
    I=mb*((l1+l2)**2.0)/12.0+mb*(Lb-l1)**2.0        #Bar MOI  
    rg=((l1**2.0-l1*l2+l2**2.0)/3)**0.5
    I2=mb*rg**2.0

    h=0.001
    newN=N+1
    
    T[0]=0.
    theta[0]=thetai
    phi[0]=pi-theta[0]
    
    hl5=l2*cos(pi-thetai)  #Part of l5 covered by the projection of l2 on l5
    if hl5-l5<0.05:
        psi[0]=arcsin(l5/l2)  #Max arm cocking angle
    else:
        psip2=arcsin((l5-hl5)/l3)  #Part of psi between horizontal and sling
        psip1=theta-pi/2       #Part of psi between horizontal and arm
        psi[0]=psip1+psip2     #Total psi


##    print 'thetai ' + str(theta[0]*180/pi)
##    print 'phii ' + str(phi[0]*180/pi)
##    print 'psii ' + str(psi[0]*180/pi)

    x3[0,0]=-l3*sin(psi[0]-theta[0])-l2*sin(theta[0])	#Cocked X-position
    x3[0,1]=-l3*cos(psi[0]-theta[0])+l2*cos(theta[0])			#Cocked Y-position
    

    x4[0,0]=l1*sin(theta[0])-l4*sin(phi[0]+theta[0])
    x4[0,1]=-l1*cos(theta[0])+l4*cos(phi[0]+theta[0])
    

    #Inital theta,phi,psi, x3, x4 check fine.




    PEcw[0]=M*g*(x4[0,1]+l5)
    PEp[0]=m*g*(x3[0,1]+l5)
    PEbar[0]=mb*g*(-(l1-l2)/2*cos(theta[0])+l5)
    
    KEcw[0]=0.0
    KEp[0]=0.0
    KEbar[0]=0.0

    KEtot[0]=KEcw[0]+KEp[0]+KEbar[0]
    PEtot[0]=PEcw[0]+PEp[0]+PEbar[0]
    
    Etot[0]=KEtot[0]+PEtot[0]


    rho=0.002378
    mu=3.7372e-7
    
    S=pi*(dia/2.0)**2.0
    flight_iter_max=1000
    maxRd=0

    
    Ten[0]=0;		
    Fb[0]=0;
    flag=0;			#Flag that indicated if projectile has lifted off ground (0=on ground)
    
    maxR=0;
    
    if runopts.gss == 0:
        print '\n'
        print '**Trebuchet Simulation**\n'
        print 'RK45 Tolerance = ' + str(Tol)
        if runopts.drag == 1:
            print 'Calculating range with drag'
        else:
            print 'Calculating range without drag'

        print ''
    flag = 0
    for i in range(1,N):

        if mod(i,100)==0 and runopts.gss==0:
            #print 'Iteration = ' + str(i) + ' Time = ' + str(T[i-1]+h)
            print 'Iteration = %(#)5.0d Sim Time = %(#2)5.3f h = %(#3)5.3g ' %{'#': i, '#2':T[i-1]+h, '#3':h}  

##            print 'h = ' + str(h)
##            print 'dthdt  ' + str(dthdt[i-1]*180/pi)
##            print 'dphidt  ' + str(dphidt[i-1]*180/pi)
##            print 'dpsidt  ' + str(dpsidt[i-1]*180/pi) 
##            print 'theta  ' + str(theta[i-1]*180/pi)
##            print 'phi  ' + str(phi[i-1]*180/pi)
##            print 'psi  ' + str(psi[i-1]*180/pi) 
##            print 'Range nd ' + str(Rnd[i-1])
##            print 'Rprime ' + str(Rprime[i-1])
##            savetxt('time.csv',T)
##            savetxt('theta.csv',theta)
##            savetxt('phi.csv',phi)
##            savetxt('psi.csv',psi)
    
        #Is projective in the air?
        
        if Fb[i-1]-m*g>0:  #Yes
            flag=1.0
        
        if i<=5:
            #print 'flag ' + str(flag)
            rkflag=1;
            rkiter=0;
            #[M,m,mb,l1,l2,l3,l4,l5,Lb,thetai,hmax,hmin,Tol,N]
            #Runge-Kutta-Fehlberg Method
            while rkflag==1.0:
                #print 'flag ' + str(flag)
                rkiter=rkiter+1;
                #print N1,N2,L,l,theta[i-1],phi[i-1],d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode
                        
                          
                k1=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j1=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f1=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                
                d1=k1/4.0
                d2=j1/4.0
                d3=f1/4.0
                
                k2=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j2=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f2=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)

                d1=3./32.*k1+9./32.*k2
                d2=3./32.*j1+9./32.*j2
                d3=3./32.*f1+9./32.*f2
                
                k3=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j3=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f3=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)

                d1=1932./2197.*k1-7200./2197.*k2+7296./2197.*k3
                d2=1932./2197.*j1-7200./2197.*j2+7296./2197.*j3
                d3=1932./2197.*f1-7200./2197.*f2+7296./2197.*f3
                
                k4=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j4=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f4=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)

                d1=439./216.*k1-8.*k2+3680./513.*k3-845./4104.*k4
                d2=439./216.*j1-8.*j2+3680./513.*j3-845./4104.*j4
                d3=439./216.*f1-8.*f2+3680./513.*f3-845./4104.*f4
                
                k5=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j5=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f5=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                
                
                d1=-8./27.*k1+2.*k2-3544./2565.*k3+1859./4104.*k4-11./40.*k5
                d2=-8./27.*j1+2.*j2-3544./2565.*j3+1859./4104.*j4-11./40.*j5
                d3=-8./27.*f1+2.*f2-3544./2565.*f3+1859./4104.*f4-11./40.*f5
                
                k6=h* dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                j6=h*dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)
                f6=h*dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1]+d1,phi[i-1]+d2,psi[i-1]+d3,psi[0],dthdt[i-1],dphidt[i-1],dpsidt[i-1],flag)          
        
                R1=1/h*abs(1./360.*k1-128./4275.*k3-2197./75240.*k4+1./50.*k5+2./55.*k6);  #Approx of global error
                R2=1/h*abs(1./360.*j1-128./4275.*j3-2197./75240.*j4+1./50.*j5+2./55.*j6);  #Approx of global error
                R3=1/h*abs(1./360.*f1-128./4275.*f3-2197./75240.*f4+1./50.*k5+2./55.*f6);  #Approx of global error
                #print 'R1 ' + str(R1)
                #print 'R2 ' + str(R2)
                #print 'R3 ' + str(R3)
                R=max(R1,R2)
                
                Rstore[i]=R;

                #print 'R ' + str(R)
##                print 'R ' + str(R) + ' TOL ' +str(Tol) + ' h ' +str(h)
                if R <= Tol or rkiter>10:
                    T[i]=T[i-1]+h
                    #print 'Time ' + str(T[i])
                    #print 'h ' + str(h)
                    #print 'RKiter ' + str(rkiter)
                    w =w +25./216.*k1+1408./2565.*k3+2197./4104.*k4-1./5.*k5
                    w2=w2+25./216.*j1+1408./2565.*j3+2197./4104.*j4-1./5.*j5
                    w3=w3+25./216.*f1+1408./2565.*f3+2197./4104.*f4-1./5.*f5
                                        
##                    print 'Rk'
##                    print 'w ' + str(w)
##                    print 'w2 ' + str(w2)
##                    print 'w3 ' + str(w3)
                    
##                    delta=0.84*(Tol/R)**(1./4.);
##                    #print 'delta ' + str(delta)
##                    if delta<=0.1:
##                        h=0.1*h;
##                    elif delta>=4.0:
##                        h=4.0*h;
##                    else:
##                        h=delta*h;
##        
##                    if h>hmax:
##                        h=hmax;
##            
##                    if h < hmin:
##                        h=hmin;
                    
                    hstore[i]=h
                    rkiterstore[i]=rkiter
                    #if rkiter>1:
                    #    print 'rkiter = ' + str(rkiter) + ' h = ' + str(h)
                    rkflag=0
##                    continue
                           

                delta=0.84*(Tol/R)**(1./4.)
                #print 'delta ' + str(delta)
                if delta<=0.1:
                    h=0.1*h;
                elif delta>=4.0:
                    h=4.*h;
                else:
                    h=delta*h;
        
                if h>hmax:
                    h=hmax;
            
                if h < hmin:
                    h=hmin

                    
##            print 'rk i=' + str(i)
##            print 'dthdt inc ' + str(w)
##            print 'dphidt inc ' + str(w2)
##            print 'dpsidt inc ' + str(w3)


        if i>4:
##                print 'in LM'
##                print h
                #Adaptive Linear Multistep Predictor Corrector scheme to get d$dts
                rkflag=1.
                rkiter=0.

                while rkflag==1:
                    rkiter=rkiter+1.0
##                    print 'h=' + str(h)
                    #Predict d$dts  with O(h^4) AB
                    dthdt[i]=  dthdt[i-1]+h*(55./24.*dthdt2[i-1]-59./24.*dthdt2[i-2]+37./24.*dthdt2[i-3]-3./8.*dthdt2[i-4])
                    dphidt[i]=dphidt[i-1]+h*(55./24.*dphidt2[i-1]-59./24.*dphidt2[i-2]+37./24.*dphidt2[i-3]-3./8.*dphidt2[i-4])
                    dpsidt[i]=dpsidt[i-1]+h*(55./24.*dpsidt2[i-1]-59./24.*dpsidt2[i-2]+37./24.*dpsidt2[i-3]-3./8.*dpsidt2[i-4])


                    
                    #Predict angles  with O(h^4) AB
##                    theta[i]=theta[i-1]+h*(55./24.*dthdt[i-1]-59./24.*dthdt[i-2]+37./24.*dthdt[i-3]-3./8.*dthdt[i-4])
##                    phi[i]=    phi[i-1]+h*(55./24.*dphidt[i-1]-59./24.*dphidt[i-2]+37./24.*dphidt[i-3]-3./8.*dphidt[i-4])
##                    psi[i]=    psi[i-1]+h*(55./24.*dpsidt[i-1]-59./24.*dpsidt[i-2]+37./24.*dpsidt[i-3]-3./8.*dpsidt[i-4])

                    #Evaluate new d$dt2s based on predicted data
                    dthdt2[i]=dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
                    dphidt2[i]=dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
                    dpsidt2[i]=dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)

                    #Correct d$dts with O(h^4) AM
                    dthdt_1=  dthdt[i-1]+h*(3./8.*dthdt2[i]+19./24.*dthdt2[i-1]-5./24.*dthdt2[i-2]+1./24.*dthdt2[i-3])  
                    dphidt_1=dphidt[i-1]+h*(3./8.*dphidt2[i]+19./24.*dphidt2[i-1]-5./24.*dphidt2[i-2]+1./24.*dphidt2[i-3])
                    dpsidt_1=dpsidt[i-1]+h*(3./8.*dpsidt2[i]+19./24.*dpsidt2[i-1]-5./24.*dpsidt2[i-2]+1./24.*dpsidt2[i-3])



##                    print 'AM3'
##                    print dthdt_1
##                    print dphidt_1
##                    print dpsidt_1
                    
                    

                    #Predict d$dts  with O(h^5) AB
                    dthdt[i]=  dthdt[i-1]+h*(1901./720.*dthdt2[i-1]-1387./360.*dthdt2[i-2]+109./30.*dthdt2[i-3]-637./360.*dthdt2[i-4]+251./720.*dthdt2[i-5])
                    dphidt[i]=dphidt[i-1]+h*(1901./720.*dphidt2[i-1]-1387./360.*dphidt2[i-2]+109./30.*dphidt2[i-3]-637./360.*dphidt2[i-4]+251./720.*dphidt2[i-5])
                    dpsidt[i]=dpsidt[i-1]+h*(1901./720.*dpsidt2[i-1]-1387./360.*dpsidt2[i-2]+109./30.*dpsidt2[i-3]-637./360.*dpsidt2[i-4]+251./720.*dpsidt2[i-5])

##                    print 'AB4'
##                    print dthdt[i]
##                    print dphidt[i]
##                    print dpsidt[i]

                    #Predict angles  with O(h^5) AB
##                    theta[i]=theta[i-1]+h*(1901./720.*dthdt[i-1]-1387./360.*dthdt[i-2]+109./30.*dthdt[i-3]-637./360.*dthdt[i-4]+251./720.*dthdt[i-5])
##                    phi[i]=   phi[i-1]+h*(1901./720.*dphidt[i-1]-1387./360.*dphidt[i-2]+109./30.*dphidt[i-3]-637./360.*dphidt[i-4]+251./720.*dphidt[i-5])
##                    psi[i]=   psi[i-1]+h*(1901./720.*dpsidt[i-1]-1387./360.*dpsidt[i-2]+109./30.*dpsidt[i-3]-637./360.*dpsidt[i-4]+251./720.*dpsidt[i-5])

                    #Evaluate new d$dt2s based on predicted data
                    dthdt2[i]=dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
                    dphidt2[i]=dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
                    dpsidt2[i]=dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i-1],phi[i-1],psi[i-1],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)


                    #Correct d$dts with O(h^5) AM
                    dthdt_2=  dthdt[i-1]+h*(251./720.*dthdt2[i]+646./720.*dthdt2[i-1]-264./720.*dthdt2[i-2]+106./720.*dthdt2[i-3]-19./720.*dthdt2[i-4])  
                    dphidt_2=dphidt[i-1]+h*(251./720.*dphidt2[i]+646./720.*dphidt2[i-1]-264./720.*dphidt2[i-2]+106./720.*dphidt2[i-3]-19./720.*dphidt2[i-4])  
                    dpsidt_2=dpsidt[i-1]+h*(251./720.*dpsidt2[i]+646./720.*dpsidt2[i-1]-264./720.*dpsidt2[i-2]+106./720.*dpsidt2[i-3]-19./720.*dpsidt2[i-4])
        


##                    print 'AM4'
##                    print dthdt_2
##                    print dphidt_2
##                    print dpsidt_2
                    
        
        
                    R1=1.0/h*abs(dthdt_2-dthdt_1)
                    R2=1.0/h*abs(dphidt_2-dphidt_1)
                    R3=1.0/h*abs(dpsidt_2-dpsidt_1)

                    

##                    print rkiter
##                    print dthdt_2
##                    print dthdt_1
##                    print dphidt_2
##                    print dphidt_1
##                    print dpsidt_2
##                    print dpsidt_1
##                
####
##                    print 'R1 ' + str(R1)
##                    print 'R2 ' + str(R2)
##                    print 'R3 ' + str(R3)
                    R=max(R1,R2)
                    Rstore[i]=R;

##
##                    print 'LM R ' + str(R) + ' TOL ' +str(Tol) + ' h ' +str(h)
                    
                    
##                    print 'LM R ' + str(R) + ' TOL ' +str(Tol) + ' h ' +str(h)
                    delta=0.84*(Tol/R)**(1./4.);
                    #print 'delta ' + str(delta)
                    if delta<=0.1:
##                        print 'Lower out'
                        h=0.1*h;
                    elif delta>=4.0:
##                        print 'up out'
                        h=4.*h;
                    else:
##                        print 'delta out'
                        h=delta*h;
            
                    #print 'New h '+ str(h)
                    if h>hmax:
                        h=hmax;
                
                    if h < hmin:
                        print 'Min h reached'
                        sys.exit()
                                    
                    
                   
                    
                                                    
                    if 1==1:  #This is where it should loop until R is < TOL but that never seems to happen.  Letting it just go works though.  Not sure why...
                        
                        

                        T[i]=T[i-1]+h
                        #Update with latest h
                        dthdt_2=dthdt[i-1]+h*(251./720.*dthdt2[i]+646./720.*dthdt2[i-1]-264./720.*dthdt2[i-2]+106./720.*dthdt2[i-3]-19./720.*dthdt2[i-4])  
                        dphidt_2=dphidt[i-1]+h*(251./720.*dphidt2[i]+646./720.*dphidt2[i-1]-264./720.*dphidt2[i-2]+106./720.*dphidt2[i-3]-19./720.*dphidt2[i-4])  
                        dpsidt_2=dpsidt[i-1]+h*(251./720.*dpsidt2[i]+646./720.*dpsidt2[i-1]-264./720.*dpsidt2[i-2]+106./720.*dpsidt2[i-3]-19./720.*dpsidt2[i-4]) 
                        #print 'Time ' + str(T[i])
                        #print 'h ' + str(h)
                        #print 'RKiter ' + str(rkiter)
                        w =dthdt_2
                        w2=dphidt_2
                        w3=dpsidt_2
        
##                        print 'dthdt inc ' + str(w)
##                        print 'dphidt inc ' + str(w2)
##                        print 'dpsidt inc ' + str(w3)
                        
##                        print 'LM'
##                        print 'w ' + str(w)
##                        print 'w2 ' + str(w2)
##                        print 'w3 ' + str(w3)

                        
                        hstore[i]=h
                        rkiterstore[i]=rkiter
                        #if rkiter>1:
                        #    print 'rkiter = ' + str(rkiter) + ' h = ' + str(h)
                        rkflag=0
##                        continue
                               








        dthdt[i]=w
        dphidt[i]=w2
       


        #Calculate new theta and phi
        if i==1:
            theta[i]=theta[i-1]+h*dthdt[i-1]  #Predictor
            phi[i]=phi[i-1]+h*dphidt[i-1]
        elif i==2:
            theta[i]=theta[i-1]+3./2.*h*dthdt[i-1]-1./2.*h*dthdt[i-2]
            phi[i]=phi[i-1]+3./2.*h*dphidt[i-1]-1./2.*h*dphidt[i-2]
        elif i==3:
            theta[i]=theta[i-1]+h*(23./12.*dthdt[i-1]-4./3.*dthdt[i-2]+5./12.*dthdt[i-3])
            phi[i]=   phi[i-1]+h*(23./12.*dphidt[i-1]-4./3.*dphidt[i-2]+5./12.*dphidt[i-3])           
        elif i==4:
            theta[i]=theta[i-1]+h*(55./24.*dthdt[i-1]-59./24.*dthdt[i-2]+37./24.*dthdt[i-3]-3./8.*dthdt[i-4])
            phi[i]=phi[i-1]+h*(55./24.*dphidt[i-1]-59./24.*dphidt[i-2]+37./24.*dphidt[i-3]-3./8.*dphidt[i-4])
        else:
            theta[i]=theta[i-1]+h*(1901./720.*dthdt[i-1]-1387./360.*dthdt[i-2]+109./30.*dthdt[i-3]-637./360.*dthdt[i-4]+251./720.*dthdt[i-5])
            phi[i]=phi[i-1]+h*(1901./720.*dphidt[i-1]-1387./360.*dphidt[i-2]+109./30.*dphidt[i-3]-637./360.*dphidt[i-4]+251./720.*dphidt[i-5])
        
        
        #print 'Phi ' + str(phi[i]*180/pi)
        


        if flag==0:
            #Since no answer for dpsidt and dpsidt2 when on ground from function
            dpsidt[i]=dthdt[i]-(l2/l3*sin(theta[i])*dthdt[i])/sqrt(1-(l2/l3*(sin(psi[0])+cos(theta[i])))**2)
            w3=dpsidt[i]
        else:
            dpsidt[i]=w3   
          
        
        if flag==0:
            psi[i]=theta[i]-pi/2+arcsin(l2/l3*(sin(psi[0])+cos(theta[i])))
        else:    
            if i==1:
                psi[i]=psi[i-1]+h*dpsidt[i-1]
            elif i==2:
                psi[i]=psi[i-1]+3./2.*h*dpsidt[i-1]-1./2.*h*dpsidt[i-2]
            elif i==3:
                psi[i]=psi[i-1]+h*(23./12.*dpsidt[i-1]-4./3.*dpsidt[i-2]+5./12.*dpsidt[i-3])
            elif i==4:
                psi[i]=psi[i-1]+h*(55./24.*dpsidt[i-1]-59./24.*dpsidt[i-2]+37./24.*dpsidt[i-3]-3./8.*dpsidt[i-4])
            else:
                psi[i]=psi[i-1]+h*(1901./720.*dpsidt[i-1]-1387./360.*dpsidt[i-2]+109./30.*dpsidt[i-3]-637./360.*dpsidt[i-4]+251./720.*dpsidt[i-5])
 
        #print 'Psi ' + str(psi[i]*180/pi)
        
        
        dthdt2[i]=dth2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i],phi[i],psi[i],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
        dphidt2[i]=dphi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i],phi[i],psi[i],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)
        dpsidt2[i]=dpsi2(M,m,mb,l1,l2,l3,l4,l5,Lb,theta[i],phi[i],psi[i],psi[0],dthdt[i],dphidt[i],dpsidt[i],flag)

    
        #Calculate projectile position
        x3[i,0]=-l3*sin(psi[i]-theta[i])-l2*sin(theta[i])
        x3[i,1]=-l3*cos(psi[i]-theta[i])+l2*cos(theta[i])
        
        #And velocity
        v3[i,0]=-l3*cos(psi[i]-theta[i])*(dpsidt[i]-dthdt[i])-l2*cos(theta[i])*dthdt[i]
        v3[i,1]= l3*sin(psi[i]-theta[i])*(dpsidt[i]-dthdt[i])-l2*sin(theta[i])*dthdt[i]
        Vtotp=sqrt(v3[i,0]**2.+v3[i,1]**2.)
        
        #And Acceleration   
        a3[i,0]=-l3*cos(psi[i]-theta[i])*(dpsidt2[i]-dthdt2[i])+(dpsidt[i]-dthdt[i])*l3*sin(psi[i]-theta[i])*(dpsidt[i]-dthdt[i])-l2*cos(theta[i])*dthdt2[i]+sin(theta[i])*dthdt[i]**2
        a3[i,1]= l3*sin(psi[i]-theta[i])*(dpsidt2[i]-dthdt2[i])+(dpsidt[i]-dthdt[i])*l3*cos(psi[i]-theta[i])*(dpsidt[i]-dthdt[i])-l2*sin(theta[i])*dthdt2[i]-cos(theta[i])*dthdt[i]**2
        a3tot=sqrt(a3[i,0]**2.+a3[i,1]**2.)

        #Calculate counter weight position
        x4[i,0]=l1*sin(theta[i])-l4*sin(phi[i]+theta[i])
        x4[i,1]=-l1*cos(theta[i])+l4*cos(phi[i]+theta[i])
        #And velocity
        v4[i,0]=l1*cos(theta[i])*dthdt[i]-l4*cos(phi[i]+theta[i])*(dphidt[i]+dthdt[i])
        v4[i,1]=l1*sin(theta[i])*dthdt[i]-l4*sin(phi[i]+theta[i])*(dphidt[i]+dthdt[i])
        VtotCW=sqrt(v4[i,0]**2.+v4[i,1]**2.)
        
        xtip[i,0]=-l2*sin(theta[i])
        xtip[i,1]=l2*cos(theta[i])
        
        x1[i,0]=l1*sin(theta[i])
        x1[i,1]=-l1*cos(theta[i])

        xbarcg[i,0]=0
        xbarcg[i,1]=0

        PEcw[i]=M*g*(x4[i,1]+l5)
        PEp[i]=m*g*(x3[i,1]+l5)
        PEbar[i]=mb*g*(-(l1-l2)/2*cos(theta[i])+l5)

        KEcw[i]=1./2.*M*VtotCW**2.
        KEp[i] =1./2.*m*Vtotp**2.
        KEbar[i]=mb/6.0*(l1**2-l1*l2+l2**2)*dthdt[i]**2

        KEtot[i]=KEcw[i]+KEp[i]+KEbar[i]
        PEtot[i]=PEcw[i]+PEp[i]+PEbar[i]
        
    
        Etot[i]=KEtot[i]+PEtot[i]



        #Calculate sling tension
        #print 'a3 ' + str(a3[i,1])
        #Ten[i]=m*a3[i,1]/sin(pi/2+psi[i]-theta[i])
        Ten[i]=m*sqrt(a3[i,0]**2+(a3[i,1])**2)
        #print 'Tension ' +str(Ten[i])

            
        #Calculate if sling has lifted projectile from ground
        Fb[i]=m*a3[i,1]
        #print 'Fb ' + str(Fb[i])

        #Calculate the direction of the projectile (in sling)

        if x3[i,0]==x3[i-1,0]:
            ang=0.0
        else:
            ang=arctan((x3[i,1]-x3[i-1,1])/(x3[i,0]-x3[i-1,0]));
            
        angout[i]=ang
        #print 'ang ' + str(ang*180/pi)
        
        #Calculate range if released from current position (assuming zero drag)
    
        Rnd[i]=(2.0*Vtotp**2.0*sin(ang)*cos(ang)/g)
        #print  'flag ' + str(flag)
        #print 'Rnd2 ' + str(Rnd[i])
        
            
        if runopts.drag==1:
            h2=0.005
        
            Vtot=sqrt(v3[i,0]**2+v3[i,1]**2)


            Re=rho*Vtot*dia/mu

        
            rvel[0,0]=v3[i,0]
            rvel[0,1]=v3[i,1]
            #print vp[i,:]   

            cd=0.27
            
            Pos[0,1]=x3[i,1]
            if rvel[0,0]>1 and rvel[0,1]>1:
                for j in range(1,flight_iter_max):
                    
                    j1=h2*accelx(rvel[j-1,0])
                    k1=h2*accely(rvel[j-1,1])

                    m1=h2*(rvel[j-1,0]+j1)
                    n1=h2*(rvel[j-1,1]+k1)


                    j2=h2*accelx(rvel[j-1,0]+0.5*j1)
                    k2=h2*accely(rvel[j-1,1]+0.5*k1)
    
                    m2=h2*(rvel[j-1,0]+0.5*j2+0.5*h2*m1)
                    n2=h2*(rvel[j-1,1]+0.5*k2+0.5*h2*n1)

                   
                    rvel[j,0]=rvel[j-1,0]+j2
                    rvel[j,1]=rvel[j-1,1]+k2

                    #print 'Vx '+ str(V[i,0])
                    #print 'Vy '+ str(V[i,1])
                   

                    Pos[j,0]=Pos[j-1,0]+m2
                    Pos[j,1]=Pos[j-1,1]+n2

         
                    Tp[j]=(j-1)*h2;
                
                    if Pos[j,1]<0:
                        #print Pos[j,1]
                        break
        
                Rwd[i]=Pos[j,0]
            else:
                Rwd[i]=0.0

        
            Rprime[i]=Rwd[i]-Rwd[i-1]

        else:  #No drag
            Rprime[i]=(Rnd[i]-Rnd[i-1])


        
        #print max(Pos[:,0])

        


        Fp[i,0]=Ten[i]*cos(psi[i])  #Direction from tip to CW end of bar is positive
        Fp[i,1]=Ten[i]*sin(psi[i])  #down is positive
        
        Farmw[i,0]=g*mb*sin(theta[i]-pi/2)
        Farmw[i,1]=g*mb*cos(theta[i]-pi/2)

        Fcw[i,0]=g*M*cos(phi[i])
        Fcw[i,1]=g*M*sin(phi[i])

        Rpin[i,0]=-(Fp[i,0]+Farmw[i,0]+Fcw[i,0])
        Rpin[i,1]=-(Fp[i,1]+Farmw[i,1]+Fcw[i,1])
        

        #Moment is incorrect !!!!!!!!!!!!!!!!!!!
        x=(l2+l1-Lb)
        Mombar1[i]=-Fp[i,1]*x+I*dthdt2[i]
        x=l2
        Mombar2[i]=-Fp[i,1]*x-Farmw[i,1]*(x+Lb-l2-l1)+I*dthdt2[i]
        x=l2+l1
        Mombar3[i]=-Fp[i,1]*x-Farmw[i,1]*(x+Lb-l2-l1)-Rpin[i,1]*(x-l2)+I*dthdt2[i]
        

            

            #Find max range to determine best release point
       

       	#Stop if arm is past set angle or sling has whirled past set angle  #theta[i]>100*pi/180 or phi[i]>270*pi/180:
        #Stop if Range is positive and derivative of Range is negative  Rnd[i]>2 and Rprime[i]<-0.005:
        if runopts.drag == 1:
            if Rwd[i]>2 and Rprime[i]<-0.1:
                newN=i
                break
        else:
            if Rnd[i]>maxR:
                maxR=Rnd[i]

            if Rnd[i]>10 and Rnd[i]<0.9*maxR:
                newN=i
                break  
    
    


##    f= open('output.csv', 'w')
##    f.write('Time th dth phi dphi psi dspi R\n')
##    for J in range(1,i):
##        f.write('%(#1)5.3f %(#2)5.3f %(#3)5.3f %(#4)5.3f %(#5)5.3f %(#6)5.3f %(#7)5.3f %(#8)5.3f %(#9)5.3f\n' %{'#1': float(T[J]),'#2': float(theta[J]*180/pi),'#3': float(dthdt[J]*180/pi),'#4': float(phi[J]*180/pi),'#5': float(dphidt[J]*180/pi),'#6': float(psi[J]*180/pi),'#7': float(dpsidt[J]*180/pi),'#8': float(Rnd[J]),'#9': float(Etot[J])})
##
##    f.close()
    

    #Find max range to determine best release point
    if runopts.drag == 1:
        [maxR,I]=maxval(Rwd)
    else:
        [maxR,I]=maxval(Rnd)
        
    if runopts.gss==0:    
        print 'Last Iteration = ' + str(I)
    
    #Find launch angle
    ang=arctan((x3[I,1]-x3[I-1,1])/(x3[I,0]-x3[I-1,0]));

    if runopts.verbose==1:
        print '\n\n***Treb Specs'
        print 'Short arm length = %(#)5.3f ft' %{'#': float(l1)}
        print 'Long arm length = %(#)5.3f ft' %{'#': float(l2)}
        print 'CW arm length = %(#)5.3f ft' %{'#': float(l4)}
        print 'Sling length = %(#)5.3f ft' %{'#': float(l3)}
        print 'Arm CG = %(#)5.3f ft' %{'#': float(Lb)}
        print 'Counter Weight = %(#)5.3f lb' %{'#': float(M*g)}
        print 'Projectile Weight = %(#)5.3f lb' %{'#': float(m*g)}
        print 'Cocked Arm Angle = %(#)5.2f deg' %{'#': float(theta[0]*180/pi)}


        print '\n\n***Sim Results'
        print 'Firing Time = %(#)5.3f s' %{'#': float(T[I])}
        
        if runopts.drag==1:
            print 'Opt Release Angle (with drag) = %(#)5.3f deg' %{'#': float(ang*180/pi)}  
            print 'Range (with drag) = %(#)5.1f ft' %{'#': float(Rwd[I])}
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])} 
        else:
            print 'Opt Release Angle (no drag) = %(#)5.3f deg' %{'#': float(ang*180/pi)}  
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])}   

        
        print 'Release Angles -- Theta = %(#)5.3f deg  -- Psi = %(#2)5.3f deg  -- Phi = %(#3)5.3f deg' %{'#': float(theta[I]*180/pi),'#2': float(psi[I]*180/pi),'#3': float(phi[I]*180/pi)}
        
        print 'Total Energy -- Initial = %(#)5.3f lb-ft  -- Max = %(#2)5.3f lb-ft   -- Min = %(#3)5.3f lb-ft   -- Delta = %(#4)5.3f lb-ft' %{'#': float(Etot[1]),'#2': float(max(Etot[0:I])),'#3': float(min(Etot[0:i])),'#4': float(max(Etot[0:I]))-float(min(Etot[0:i]))}
        print 'Energy Efficiency into projectile = %(#)5.3f %%' %{'#': float(KEp[I]/(Etot[0]-PEcw[I])*100)}



        elapsed = time.time()-start
        print 'Program run time = %(#)5.3f s' %{'#': float(elapsed)}
    
    if runopts.verbose==0 and runopts.gss==0:
        if runopts.drag==1:
            print 'Range (with drag) = %(#)5.1f ft' %{'#': float(Rwd[I])}
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])} 
        else:
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])}  
    
    if runopts.interactive==1:
        opt=1
        
        while (opt != 0):
            print '\n\n'
            print ' 1 - Range'
            print ' 2 - Angles'
            print ' 3 - Sling Tension'
            print ' 4 - Energy'
            print ' 5 - h values'
            print ' 6 - R values'
            print ' 7 - Arm Moment'
            print ' 8 - Strobe of Firing'
            print ' 9 - Derivative Values'
            print ' 0 - Quit'

            
            try:
                opt = int(raw_input('Selection: '))
            except ValueError:
                print '***Invalid choice***'
                continue
            
            
            if opt == 1:
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                if runopts.drag==1:
                    ax1.plot(T[0:i],Rnd[0:i],label='Range',color='blue')
                    ax1.plot(T[0:i],Rwd[0:i],label='Range wd',color='green')
                    ax1.plot([T[I], T[I]],[Rwd[I], Rwd[I]],'x',label='Launch Point',color='red',mew=2)
                else:
                    ax1.plot(T[0:i],Rnd[0:i],label='Range',color='blue')
                    ax1.plot([T[I], T[I]],[Rnd[I], Rnd[I]],'x',label='Launch Point',color='red',mew=2)
                ax1.set_xlabel('Time (s)')
                ax1.set_ylabel('Range (ft)')
                    
                ax2 = ax1.twinx()    
                ax2.plot(T[0:i],Rprime[0:i],label='Rprime',color='black')
                
                ax1.grid(True)
                ax1.legend(loc=2)
                ax2.legend(loc=4)
                plt.show()
            elif opt == 2:
                plt.plot(T[0:I+1],theta[0:I+1]*180/pi,label='Theta')
                plt.plot(T[0:I+1],phi[0:I+1]*180/pi,label='Phi')
                plt.plot(T[0:I+1],psi[0:I+1]*180/pi,label='Psi')
                plt.plot(T[0:I+1],angout[0:I+1]*180/pi,label='Proj Angle')
                plt.xlabel('Time (s)')
                plt.ylabel('Angle (deg)')
                plt.grid(True)
                plt.legend(loc=2)
                plt.show()
            elif opt == 3:
                plt.plot(T[0:I+1],Ten[0:I+1],color='blue')
                plt.xlabel('Time (s)')
                plt.ylabel('Tension (lbf)')
                plt.grid(True)
                plt.show()   
            elif opt == 4:
                plt.plot(T[0:I+1],PEcw[0:I+1],label='PEcw')   
                plt.plot(T[0:I+1],PEp[0:I+1],label='PEp')
                plt.plot(T[0:I+1],PEbar[0:I+1],label='PEbar')
                plt.plot(T[0:I+1],KEcw[0:I+1],label='KEcw') 
                plt.plot(T[0:I+1],KEp[0:I+1],label='KEp')
                plt.plot(T[0:I+1],KEbar[0:I+1],label='KEbar') 
                plt.plot(T[0:I+1],Etot[0:I+1],label='Etot')
                plt.xlabel('Time (s)')
                plt.ylabel('Energy')
                plt.grid(True)
                plt.legend(loc=6)
                plt.show()
            elif opt==5:
                plt.plot(hstore[0:I+1])
                plt.xlabel('Iteration')
                plt.ylabel('h')
                plt.grid(True)
                plt.show()
            elif opt==6:
                plt.plot(Rstore[0:I+1])
                plt.xlabel('Iteration')
                plt.ylabel('R')
                plt.grid(True)
                plt.show()
            elif opt==7:
                #plt.plot(T[0:I+1],Rpin[0:I+1,0], color='blue')
                #plt.plot(T[0:I+1],Rpin[0:I+1,1], color='red')
                #plt.plot(T[0:I+1],Mombar1[0:I+1], color='blue')
                #plt.plot(T[0:I+1],Mombar2[0:I+1], color='red')
                #plt.plot(T[0:I+1],Mombar3[0:I+1], color='green')
                #plt.plot(Fyp3[0:I+1])
                plt.plot(T[0:I+1],a3[0:I+1,1])
                plt.xlabel('Arm location')
                plt.ylabel('Moment')
                plt.grid(True)
                plt.show()
            elif opt==8:
                j=floor(linspace(1,I+1,10, endpoint=True))
                for i in j:
                    plt.plot([x1[i,0], xtip[i,0]],[x1[i,1], xtip[i,1]], color='blue')
                    plt.plot([x3[i,0], xtip[i,0]],[x3[i,1], xtip[i,1]], color='red')
                    plt.plot([x1[i,0], x4[i,0]],[x1[i,1], x4[i,1]], color='red')
                    plt.plot([x4[i,0],x4[i,0]],[x4[i,1],x4[i,1]],marker='o',color='green')
                    #plt.plot([-l*cos(theta[i]), -l*cos(theta[i])],[0,0],marker='o',color='yellow')
                    
                plt.grid(True)
                plt.show() 
            elif opt==9:
                plt.plot(T[0:I+1],dthdt[0:I+1],label='dthdt',color='blue')
                plt.plot(T[0:I+1],dphidt[0:I+1],label='dphidt',color='green')
                plt.plot(T[0:I+1],dpsidt[0:I+1],label='dpsidt',color='red')
                plt.xlabel('Time')
                plt.ylabel('rad/s')
                plt.legend(loc=2)
                plt.grid(True)
                plt.show()
            elif opt==0:
                sys.exit()

    
    else:
        if runopts.drag==1:
            return Rwd[I]
        else:
            return Rnd[I]
          
                




def func(inputs,x1,s,AL):
    #print '\nx1 = ' + str(x1)

    #print s
    #print x1
    #[M,m,mb,l1,l2,l3,l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]=inputs
    inputs[3]=float(inputs[3]+x1*s[0])
    inputs[4]=float(AL-inputs[3])
    inputs[5]=float(inputs[5]+x1*s[1])    
    
    #print 'Trying l1 = '+ str(inputs[2]) + ' l2 = ' + str(inputs[3]) + ' d = ' + str(inputs[4])
    if inputs[2]<0.1 or inputs[3]<0.1 or inputs[4]<0:
        out = -(abs(inputs[2])+abs(inputs[3])+abs(inputs[4]))**2
    else:
        out =  main(inputs,runopts)
    #print '*Func Range = ' + str(out)
   
    return out


def gss(inputs,a,b,eps,N,s,AL):
    #GSS to maximize range
    C=(-1+sqrt(5))/2
    x1=C*a+(1-C)*b
    print 'fx1'
    fx1 = float(func(inputs[:],x1,s,AL))
    x2 = (1-C)*a + C*b
    print 'fx2'
    fx2 = float(func(inputs[:],x2,s,AL))
    print '-----------------------------------------------------------------------------------------------------------------'
    print '       a             x1                x2                b          f(x1)           f(x2)           b - a'
    print '-----------------------------------------------------------------------------------------------------------------'
    print ' %(#)9.5g\t %(#2)9.5g\t %(#3)9.9g\t %(#4)9.9g\t  %(#5)9.9g   %(#6)9.9g   %(#7)9.9g' %{'#': float(a),'#2':float(x1),'#3':float(x2),'#4':float(b),'#5':float(fx1),'#6':float(fx2),'#7':float(b-a)}
    #print str(x1)+'\t'+str(x2)+'\t'+str(fx1)+'\t'+str(fx2)+'\t'+str(x2-x1)
    
    for i in range(0,N):
        #print 'fx1 ' + str(fx1) + ' fx2 ' + str(fx2)
        if fx1 > fx2:
            #print 'fx1 larger'
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = C*a + (1-C)*b
            fx1 = float(func(inputs[:],x1,s,AL))
        else:
            #print 'fx2 larger'
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = (1-C)*a + C*b
            fx2 = float(func(inputs[:],x2,s,AL))        
         
        print ' %(#)9.5g\t %(#2)9.5g\t %(#3)9.9g\t %(#4)9.9g\t  %(#5)9.9g   %(#6)9.9g   %(#7)9.9g' %{'#': float(a),'#2':float(x1),'#3':float(x2),'#4':float(b),'#5':float(fx1),'#6':float(fx2),'#7':float(b-a)}   
        if (abs(b-a) < eps) and abs(fx1-fx2)<1:
            print('succeeded after '+str(i)+ ' steps');
            out=(b+a)/2;
            return out


    print('failed requirements after max steps');
        
    


                



















#cProfile.run('main()',sort=1, filename="treb_adapt.cprof")

#stats = pstats.Stats("treb_adapt.cprof")
#stats.print_stats()
#stats.strip_dirs().sort_stats('time').print_stats(20)





















if __name__ == "__main__":
    startt=time.time()
    runopts = opts(0,1)
    inputs = getinputfromfile()
    step=0.01
    endtol=0.5
    x=array([0.0,0.0])
    z=array([0.0,0.0])
    A=mat(zeros([3,3]))
    B=mat(zeros([3,1]))
    #AL=3.0 				
    xstore=zeros([20,2])
    astore=zeros([20,3])
    stopflag=0
    
    for i in range(0,20):  #while resid>0.1:
        print '\n\n********* STARTING ITERATION ' + str(i+1) + ' ***************\n\n'
        [M,m,mb,l1,l2,l3,l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]=inputs

        AL=l1+l2  #Set the arm length to the sum of the parts
        
        if i>0:
            base=gtnew
        else:
            print '\n**Calculating intial range.'
            base=main(inputs,runopts)
            
        x[0]=l1+step  #vary short part of the arm (axle to tip)
        x[1]=l3+step   #vary sling length
        
        inputs=[M,m,mb,x[0],AL-x[0],l3,l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]  #range with varied arm ratio
        print '\n**Calculating range with change in arm ratio.'
        z[0]=main(inputs,runopts)
        
        inputs=[M,m,mb,l1,l2,x[1],l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]  #range with varied sling length
        print '\n**Calculating range with change in sling length.'
        z[1]=main(inputs,runopts)
        print 'z = ' + str(z-base)

        z=(z-base)/step  #gradient


        
        z0=sqrt(z[0]**2+z[1]**2)
        print 'z0 = ' + str(z0)
        z=z/z0   #Normalize gradient
        print 'norm z = ' + str(z)
        #time.sleep(3)
        
        print 'In GSS'
        runopts.gss=1
        alpha=gss(inputs,0,1,10**-3.,50,z,AL)
        runopts.gss=0
    
        
        print 'alpha = ' + str(alpha)
        
        x=array(x+alpha*z)

        print 'xchange = ' + str(+alpha*z)
        print 'New x ' + str(x)
        
        #time.sleep(10)

        inputs=[M,m,mb,x[0],AL-x[0],x[1],l4,l5,Lb,thetai,dia,hmax,hmin,Tol,N]
        print 'Calculating new range.'
        if runopts.verbose == 0:
            runopts.verbose = 1
            gtnew=main(inputs,runopts)
            runopts.verbose = 0
        else:
            gtnew=main(inputs,runopts)


        print 'base range = ' + str(base)
        print 'gtnew range = ' + str(gtnew) + ' alpha = ' + str(alpha)
        
        resid=gtnew-base
        if abs(resid)<=1.0:
                print 'Optimum found!'
                stopflag=1
                
        
        
        print 'Increase in range = ' + str(resid) + ' ft'
        print '****New Range = ' + str(gtnew) + ' ft'
        if abs(resid)<endtol:
            Imax=i

            elapsed = time.time()-startt
            print 'Program run time = %(#)5.3f min' %{'#': float(elapsed/60)}
            break

        if stopflag==1:
            sys.exit()

        #time.sleep(5)
        
    #plt.plot(xstore[0:Imax,0],xstore[0:Imax,1])
    #plt.show()

    #plt.plot(astore[0:Imax,0],color='blue')
    #plt.plot(astore[0:Imax,1],color='green')
    #plt.plot(astore[0:Imax,2],color='red')
    #plt.show()
        
    
