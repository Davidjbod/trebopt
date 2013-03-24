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
    def __init__(self, opt1, opt2, opt3):
        self.interactive = opt1
        self.drag = opt2
        self.verbose = opt3



def dth2(N1,N2,L,l,theta,phi,d,g,M,m,mb,I,Lb,dthdt,dphidt,mode):

    if mode == 0:
    	A=1/N1*m*((L+l)*sin(theta)*cos(phi)-L*cos(theta)*sin(phi));
    	B=-((L+l)*cos(theta)*cos(phi)+L*sin(theta)*sin(phi))*dthdt**2+d*dphidt**2+N2;
    	C=M*l*cos(theta)*(l*sin(theta)*dthdt**2+g);
    	D=(M*l**2*cos(theta)**2+(m*((L+l)*sin(theta)*cos(phi)-L*cos(theta)*sin(phi))**2)/N1)**-1;

    	out=(A*B+C)*D;
    else:	
        A=(L+l)*sin(theta)*cos(phi)-L*cos(theta)*sin(phi)

        A1=cos(theta)*(A*m*(L+l)+l*sin(theta)*cos(phi)*(-M*l+mb*(Lb-l)+mb*Lb))
        B1=sin(theta)*sin(phi)/cos(phi)*(l*cos(theta)*sin(phi)*(-M*l+mb*(Lb-l)+mb*Lb)+L*A*m)
        C1=dphidt**2*A*m*d*(cos(phi)+sin(phi)**2/cos(phi))
        D1=g*cos(theta)*(mb*Lb-M*l-mb*l)*(sin(phi)**2/cos(phi)+cos(phi))
        E1=-cos(phi)*(M*l**2*cos(theta)**2-mb*(Lb-l)*l*cos(theta)**2+mb*Lb*l*sin(theta)**2+I)-A*m*(L+l)*sin(theta)
        F1=sin(phi)/cos(phi)*(L*A*m*cos(theta)-sin(phi)*(M*l**2*cos(theta)**2-mb*(Lb-l)*l*cos(theta)**2+mb*Lb*l*sin(theta)**2+I))

        out=(dthdt**2*(A1+B1)-C1+D1-A*m*g*sin(phi)/cos(phi))/(E1+F1)
    return out


def dphi2(N1,N2,L,l,theta,phi,d,g,M,m,mb,I,Lb,dthdt,dphidt,dth2dt2,flag,mode):
    
    if flag==1:
        if mode == 0:
            A=((L+l)*sin(theta)*sin(phi)+L*cos(theta)*cos(phi))*dth2dt2
            B=((L+l)*cos(theta)*sin(phi)-L*sin(theta)*cos(phi))*dthdt**2+g*cos(phi)
            out=1/d*(A+B)
        else:
            A=(L+l)*sin(theta)*cos(phi)-L*cos(theta)*sin(phi)
            
            A1=dth2dt2*(L*A*m*cos(theta)-sin(phi)*(M*l**2*cos(theta)**2-mb*(Lb-l)*l*cos(theta)**2+mb*Lb*l*sin(theta)**2+I))
            B1=-dthdt**2*sin(theta)*(l*cos(theta)*sin(phi)*(-M*l+mb*(Lb-l)+mb*Lb)+L*A*m)
            C1=A*d*m*sin(phi)*dphidt**2
            D1=-g*sin(phi)*cos(theta)*(mb*Lb-M*l-mb*l)+A*m*g

            out=(A1+B1+C1+D1)/(A*m*d*cos(phi))
    else:
        out = 1/(d*cos(phi))*(L*cos(theta)*dth2dt2-L*sin(theta)*dthdt**2+d*sin(phi)*dphidt**2);   

    return out



#Find max val and location of max val
def maxval(arr):
    val=-10**12
    index=0
    for i in range(0,len(arr)):
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
def getinputfromfile(runopts):
    

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
    M=float(indata[5])/g;      #Counter weight mass - slugs
    m=float(indata[6])/g;        #Projectile mass - slugs

    L=float(indata[1])	      #Distance from tip of throwing arm to wheel axle - ft
    l=float(indata[2])	      #Distance from counter weight cg to wheel axle - ft
    d=float(indata[3])            #Sling Length - ft
    th=float(indata[4])              #Track Height - ft

    mb=float(indata[7])/g         #Bar mass
    if float(indata[8])==-1.0:
        Lb=(L+l)/3.         #Bar CG
    else:
        Lb=indata[8]

    #Sim Setup
    hmax=float(indata[14])       #Max dt for Runge-Kutta Method
    hmin=float(indata[15])       #Min dt for RK
    Tol=float(indata[16])        #Tolerance for R error
    N=int(indata[17])  	      #Max number of iterations

    runopts.interactive = int(indata[19])
    runopts.drag = int(indata[20])
    runopts.verbose = int(indata[21])
                              

    #Inital values
    thetai=float(indata[9])*pi/180.	
    phii=float(indata[10])*pi/180.

    dia=float(indata[11])
    mode=int(indata[13]);                     #Original equation = 0,  My Equations = 1

    return [M,m,L,l,d,th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]

#Main code
def main(inputs,runopts):
    global rho
    global g
    global m
    global S
    global cd
    
    start = time.time()
    

    [M,m,L,l,d,th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]=inputs

    if L<=0.1 or l<=0.1 or d<=0.1:
        return 0
    

    w=0               #Internal variable for RK
    w2=0              #Internal variable for RK

    #Initialize arrays and matricies

    xp=zeros([N,2])	#projectile position in sling
    vp=zeros([N,2])	#projective velocity in sling
    xcw=zeros([N,2])	#Counter weight position
    vcw=zeros([N,2])	#Counter weight velocity
    xtip=zeros([N,2])
    xbarcg=zeros([N,2])
    
    dthdt=zeros([N,1])	#Derivative of Theta wrt time
    dphidt=zeros([N,1])	#Derivative of Phi wrt time
    theta=zeros([N,1])	#Arm Angle (0 deg = horizonal, Neg. angles = CW above track
    phi=zeros([N,1])	#Sling Angle (0 deg = horizontal, Pos. angles = normal movement *180/pi
    T=zeros([N,1])	#Time - s
    Ten=zeros([N,1])	#Sling Tension - lbf
    Fb=zeros([N,1])	#Force on projectile in the vertical direction (Pos. = up) - lbf
    Rnd=zeros([N,1])	#Calculated zero drag range if released - ft

    PEcw=zeros([N,1])	#Potential Energy of Counter weight
    PEp =zeros([N,1])	#Potential Energy of Projectile
    KEcw=zeros([N,1])	#Kenetic Energy of Counter weight
    KEp=zeros([N,1])	#Kenetic Energy of Projectile
    PEbar=zeros([N,1])	#Potenial Energy of bar
    KEbar=zeros([N,1])	#Kenetic Energy of bar    
    Etot=zeros([N,1])       #Total Energy

    Rstore=zeros([N,1])     #RK Error diff
    hstore=zeros([N,1])     #h storage
    rkiterstore=zeros([N,1])     #h storage

    Fyp1=zeros([N,1])       #Bar force in yp due to M*g
    Fyp2=zeros([N,1])       #Bar force in yp due to RM
    Fyp3=zeros([N,1])       #Bar force in yp due to Ra
    Fyp4=zeros([N,1])       #Bar force in yp due to Mb*g
    Fyp5=zeros([N,1])       #Bar force in yp due to T

    Mombar13=zeros([N,1])
    Mombar34=zeros([N,1])
    Mombar45=zeros([N,1])

    angout=zeros([N,1])
    Rprime=zeros([N,1])
    Vtot=zeros([N,1])

    out=zeros([1000,2])
    rvel=zeros([1000,2])
    Pos=zeros([1000,2])
    Rwd=zeros([N,1])
    Tp=zeros([1000,1])
    
  
    


    g=32.17;      #Gravity - ft/s^2
    I=mb*(Lb-l)**2.+M*l**2.          #Bar MOI
    
    h=hmax
   	
    newN=N+1
		
    T[0]=0.
    theta[0]=thetai
    phi[0]=phii

    xp[0,0]=-l*cos(theta[0])-L*cos(theta[0])+d*cos(phi[0]);	#Cocked X-position
    xp[0,1]=L*sin(theta[0])-d*sin(phi[0]);			#Cocked Y-position


    xcw[0,1]=-l*sin(theta[0])
    PEcw[0]=M*g*(xcw[0,1]-(-l))
    PEp[0]=0
    PEbar[0]=0

    KEcw[0]=0.0
    KEp[0]=0.0
    KEbar[0]=0.0

    
    Etot[0]=PEcw[0]





    rho=0.002378
    mu=3.7372e-7
    
    S=pi*(dia/2.0)**2.0
    flight_iter_max=1000
    maxRd=0

    
    Ten[0]=0;		
    Fb[0]=0;
    flag=0;			#Flag that indicated if projectile has lifted off ground (0=on ground)
    

    if runopts.verbose==1:
        print '\n\n\n'
        print '**Floating Arm Treb Simulation**\n'
        print 'RK45 Tolerance = ' + str(Tol)


        
    if runopts.drag == 1:
        if runopts.verbose==1:
            print 'Calculating range with drag'
    else:
        if runopts.verbose==1:
            print 'Calculating range without drag'

    
    for i in range(1,N):

        if mod(i,100)==0 and runopts.verbose==1:
            #print 'Iteration = ' + str(i) + ' Time = ' + str(T[i-1]+h)
            print 'Iteration = %(#)5.0d Sim Time = %(#2)5.3f' %{'#': i, '#2':T[i-1]+h}
    
        #Is projective in the air?
        if Fb[i-1]>0 or flag==1:
            #Yes
            N1=1;
            N2=g*sin(phi[i-1]);
            flag=1;
        else:
            #No
            N1=cos(phi[i-1])**2;
            N2=0;
            flag=0;
    
    
        rkflag=1;
        rkiter=0;
    	#Runge-Kutta-Fehlberg Method
        while rkflag==1:
            rkiter=rkiter+1;
            #print N1,N2,L,l,theta[i-1],phi[i-1],d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode
            dth2dt2=dth2(N1,N2,L,l,theta[i-1],phi[i-1],d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k1=h*dth2dt2
            j1=h*dphi2(N1,N2,L,l,theta[i-1],phi[i-1],d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)

            dth2dt2=dth2(N1,N2,L,l,theta[i-1]+k1/4.,phi[i-1]+j1/4.,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k2=h*dth2dt2
            j2=h*dphi2(N1,N2,L,l,theta[i-1]+k1/4.,phi[i-1]+j1/4.,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)

            dth2dt2=dth2(N1,N2,L,l,theta[i-1]+3./32.*k1+9./32.*k2,phi[i-1]+3./32.*j1+9./32.*j2,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k3=h*dth2dt2
            j3=h*dphi2(N1,N2,L,l,theta[i-1]+3./32.*k1+9./32.*k2,phi[i-1]+3./32.*j1+9./32.*j2,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)

            dth2dt2=dth2(N1,N2,L,l,theta[i-1]+1932./2197.*k1-7200./2197.*k2+7296./2197.*k3,phi[i-1]+1932./2197.*j1-7200./2197.*j2+7296./2197.*j3,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k4=h*dth2dt2
            j4=h*dphi2(N1,N2,L,l,theta[i-1]+1932./2197.*k1-7200./2197.*k2+7296./2197.*k3,phi[i-1]+1932./2197.*j1-7200./2197.*j2+7296./2197.*j3,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)

            dth2dt2=dth2(N1,N2,L,l,theta[i-1]+439./216.*k1-8.*k2+3680./513.*k3-845./4104.*k4,phi[i-1]+439./216.*j1-8.*j2+3680./513.*j3-845./4104.*j4,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k5=h*dth2dt2
            j5=h*dphi2(N1,N2,L,l,theta[i-1]+439./216.*k1-8.*k2+3680./513.*k3-845./4104.*k4,phi[i-1]+439./216.*j1-8.*j2+3680./513.*j3-845./4104.*j4,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)

            dth2dt2=dth2(N1,N2,L,l,theta[i-1]-8./27.*k1+2.*k2-3544./2565.*k3+1859./4104.*k4-11./40.*k5,phi[i-1]-8./27.*j1+2.*j2-3544./2565.*j3+1859./4104.*j4-11./40.*j5,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],mode)
            k6=h*dth2dt2
            j6=h*dphi2(N1,N2,L,l,theta[i-1]-8./27.*k1+2.*k2-3544./2565.*k3+1859./4104.*k4-11./40.*k5,phi[i-1]-8./27.*j1+2.*j2-3544./2565.*j3+1859./4104.*j4-11./40.*j5,d,g,M,m,mb,I,Lb,dthdt[i-1],dphidt[i-1],dth2dt2,flag,mode)        

    
            R=1/h*abs(1./360.*k1-128./4275.*k3-2197./75240.*k4+1./50.*k5+2./55.*k6);


            Rstore[i]=R;

        
            if R <= Tol or rkiter>10:
                T[i]=T[i-1]+h
                w=w+25./216.*k1+1408./2565.*k3+2197./4104.*k4-1./5.*k5
                w2=w2+25./216.*j1+1408./2565.*j3+2197./4104.*j4-1./5.*j5
                hstore[i]=h
                rkiterstore[i]=rkiter
                #if rkiter>1:
                #    print 'rkiter = ' + str(rkiter) + ' h = ' + str(h)
                rkflag=0
                continue
                       
            delta=0.84*(Tol/R)**(1./4.);

            if delta<=0.1:
              	h=0.1*h;
            elif delta>=4:
               	h=4.*h;
            else:
               	h=delta*h;
    
            if h>hmax:
       	       	h=hmax;
        
            if h < hmin:
                h=hmin;

   
 
        dthdt[i]=w
        dphidt[i]=w2

    
        #Calculate new theta and phi
        theta[i]=theta[i-1]+3./2.*h*dthdt[i]-1./2.*h*dthdt[i-1];
        phi[i]=phi[i-1]+3./2.*h*dphidt[i]-1./2.*h*dphidt[i-1];
    

        #Calculate projectile position
        xp[i,0]=-(L+l)*cos(theta[i])+d*cos(phi[i]);
        xp[i,1]=L*sin(theta[i])-d*sin(phi[i]);
      	#And velocity
        vp[i,0]=(L+l)*sin(theta[i])*dthdt[i]-d*sin(phi[i])*dphidt[i]
        vp[i,1]=L*cos(theta[i])*dthdt[i]-d*cos(phi[i])*dphidt[i]
        Vtot[i]=sqrt(vp[i,0]**2.+vp[i,1]**2.)

        
        #Calculate counter weight position
        xcw[i,0]=0
        xcw[i,1]=-l*sin(theta[i])
      	#And velocity
        vcw[i,0]=0
        vcw[i,1]=-l*cos(theta[i])*dthdt[i]

        xtip[i,0]=-(L+l)*cos(theta[i])
        xtip[i,1]=L*sin(theta[i])

        xbarcg[i,0]=-Lb*cos(theta[i])
        xbarcg[i,1]=(Lb-l)*sin(theta[i])

        PEcw[i]=M*g*(xcw[i,1]-(-l))
        PEp[i]=m*g*(xp[i,1]-xp[0,1])
        PEbar[i]=mb*g*((Lb-l)*sin(theta[i])-(Lb-l)*sin(theta[0]))

        KEcw[i]=1./2.*M*vcw[i,1]**2.
        KEp[i] =1./2.*m*Vtot[i]**2.
        vtotbar=sqrt((Lb*sin(theta[i])*dthdt[i])**2.+((Lb-l)*cos(theta[i])*dthdt[i])**2.)
        KEbar[i]=1./2.*mb*vtotbar**2.+1./2.*I*dthdt[i]**2.

        Etot[i]=PEcw[i]+PEp[i]+KEcw[i]+KEp[i]+PEbar[i]+KEbar[i]


        #Calculate d_theta^2/dt^2 for Tension equation
        dth2dten=dth2(N1,N2,L,l,theta[i],phi[i],d,g,M,m,mb,I,Lb,dthdt[i],dphidt[i],mode);


        #Calculate sling tension


        if mode==0:
            Ten[i]=((M*l**2*cos(theta[i])**2*dth2dten-M*l*cos(theta[i])*(l*sin(theta[i])*dthdt[i]**2+g))/((L+l)*sin(theta[i])*cos(phi[i])-L*cos(theta[i])*sin(phi[i])))
        else:
            A2=dth2dten*(M*l**2.*cos(theta[i])**2.-mb*(Lb-l)*l*cos(theta[i])**2.+mb*Lb*l*sin(theta[i])**2.+I)
            B2=dthdt[i]**2.*l*cos(theta[i])*sin(theta[i])*(-M*l+mb*(Lb-l)+mb*Lb)
            C2=g*cos(theta[i])*(mb*Lb-M*l-mb*l)
            D2=(L+l)*sin(theta[i])*cos(phi[i])-L*cos(theta[i])*sin(phi[i])
            Ten[i]=(A2+B2+C2)/D2
            
        #Calculate if sling has lifted projectile from ground
        Fb[i]=(Ten[i]*sin(phi[i])-m*g)

        #Calculate the direction of the projectile (in sling)
        ang=arctan((xp[i,1]-xp[i-1,1])/(xp[i,0]-xp[i-1,0]));
        angout[i]=ang
        #Calculate range if released from current position (assuming zero drag)
    
        Rnd[i]=(2.*Vtot[i]**2.*sin(ang)*cos(ang)/g)
        
        if Rnd[i]>1:   #Add extra range due to being released about ground level
            hstart=th+L*sin(theta[i])+d*sin(phi[I]-pi)  #Height above ground that proj. is released
            Rnd[i]=Rnd[i]+hstart/tan(ang)  #Impact angle same as release angle (approx as right angle triangle)

            
        if runopts.drag==1:
            h2=0.005
        
            Vtotdrag=sqrt(vp[i,0]**2+vp[i,1]**2)


            Re=rho*Vtotdrag*dia/mu

        
            rvel[0,0]=vp[i,0]
            rvel[0,1]=vp[i,1]
            #print vp[i,:]   

            cd=0.27
            
            Pos[0,1]=th+L*sin(theta[i])+d*sin(phi[I]-pi)
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
                        J=j
                        break
        
                Rwd[i]=Pos[j,0]
            else:
                Rwd[i]=0.0

        
            Rprime[i]=Rwd[i]-Rwd[i-1]

        else:  #No drag
            Rprime[i]=(Rnd[i]-Rnd[i-1])


        
        #print max(Pos[:,0])

        

        Fyp1[i]=-M*g*cos(theta[i])
        
        RM=-Ten[i]*cos(phi[i])
        Fyp2[i]=RM*sin(theta[i])
        
        d2yMdt=(Lb-l)*(cos(theta[i])*dth2dten-dthdt[i]*sin(theta[i]))
        Ra=M*d2yMdt+M*g+Ten[i]*sin(phi[i])
        Fyp3[i]=Ra*cos(theta[i])

        Fyp4[i]=-mb*g*cos(theta[i])
        
        Fyp5[i]=-Ten[i]*cos(pi/2.-(phi[i]-theta[i]))

        Mombar13[i]=((Fyp2[i]+Fyp1[i])-Fyp3[i])/l
        Mombar34[i]=(Fyp3[i]-Fyp4[i])/(Lb-l)
        Mombar45[i]=(Fyp4[i]-Fyp5[i])/(L-Lb)



        

       	#Stop if arm is past set angle or sling has whirled past set angle  #theta[i]>100*pi/180 or phi[i]>270*pi/180:
        #Stop if Range is positive and derivative of Range is negative  Rnd[i]>2 and Rprime[i]<-0.005:
        if runopts.drag == 1:
            if Rwd[i]>2 and Rprime[i]<-0.5:
                newN=i
                break
        else:
            if Rnd[i]>2 and Rprime[i]<-0.5:
                newN=i
                break  
    

    #Find max range to determine best release point
    if runopts.drag == 1:
        [maxR,I]=maxval(Rwd)
    else:
        [maxR,I]=maxval(Rnd)
    


    if runopts.verbose==1:
        print 'Last Iteration = ' + str(I)
    
        #Find launch angle
        ang=arctan((xp[I,1]-xp[I-1,1])/(xp[I,0]-xp[I-1,0]));

        print '\n\n***Treb Specs'
        print 'Counter weight to axle distance = %(#)5.3f ft' %{'#': float(l)}
        print 'Arm Tip to axle distance = %(#)5.3f ft' %{'#': float(L)}
        print 'Sling distance = %(#)5.3f ft' %{'#': float(d)}
        print 'Arm CG = %(#)5.3f ft' %{'#': float(Lb)}
        print 'Counter Weight = %(#)5.3f lb' %{'#': float(M*g)}
        print 'Projectile Weight = %(#)5.3f lb' %{'#': float(m*g)}
        print 'Cocked Arm Angle = %(#)5.2f deg' %{'#': float(theta[0]*180/pi)}


        print '\n\n***Sim Results'
        print 'Firing Time = %(#)5.3f s' %{'#': float(T[I])}

        if runopts.drag==1:
            print 'Opt Release Angle (with drag) = %(#)5.3f deg' %{'#': float(ang*180/pi)}  
            print 'Range (with drag) = %(#)5.5f ft' %{'#': float(Rwd[I])}
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])} 
        else:
            print 'Opt Release Angle (no drag) = %(#)5.3f deg' %{'#': float(ang*180/pi)}  
            print 'Range (no drag) = %(#)5.1f ft' %{'#': float(Rnd[I])}   


        print 'Release Angles -- Theta = %(#)5.3f deg  -- Phi = %(#2)5.3f deg' %{'#': float(theta[I]*180/pi),'#2': float(phi[I]*180/pi)}

        print 'Total Energy -- Initial = %(#)5.3f lb-ft  -- Max = %(#2)5.3f lb-ft   -- Min = %(#3)5.3f lb-ft' %{'#': float(Etot[1]),'#2': float(max(Etot[0:I])),'#3': float(min(Etot[0:i]))}
        print 'Energy Efficiency = %(#)5.3f %%' %{'#': float(KEp[I]/PEcw[1]*100)}
        Emaxmindiff=max(max(Etot[0:I])-Etot[1],Etot[1]-min(Etot[1:I]))
        print 'E Tot diff = %(#)5.3f %%' %{'#': float((Emaxmindiff)/Etot[1]*100)}


        elapsed = time.time()-start
        print 'Program run time = %(#)5.3f s' %{'#': float(elapsed)}



        if runopts.interactive==1:
            opt=1

            while (opt != 0):
                print '\n\n'
                print ' 1 - Range'
                print ' 2 - Angles'
                print ' 3 - Sling Tension'
                print ' 4 - Energy'
                print ' 5 - h values'
                print ' 6 - rkiter values'
                print ' 7 - Arm Moment'
                print ' 8 - Strobe of Firing'
                print ' 9 - Projectile Velocity'
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
                    plt.plot(T[0:I+1],angout[0:I+1]*180/pi,label='Proj Angle')
                    plt.xlabel('Time (s)')
                    plt.ylabel('Angle (deg)')
                    plt.grid(True)
                    plt.legend(loc=2)
                    plt.show()
                elif opt == 3:
                    plt.plot(T[0:I+1],Ten[0:I+1])
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
                    plt.plot(rkiterstore[0:I+1])
                    plt.xlabel('Iteration')
                    plt.ylabel('rkiter')
                    plt.grid(True)
                    plt.show()
                elif opt==7:
                    minv=min(Mombar13[0:I+1])
                    plt.plot([0, -l],[minv, minv], color='blue')
                    maxv=max(Mombar34[0:I+1])
                    plt.plot([-l, -Lb],[maxv, maxv], color='blue')
                    plt.plot([-l, -l],[minv, maxv], color='blue')
                    maxv2=max(Mombar45[0:I+1])
                    plt.plot([-Lb, -Lb],[maxv, maxv2], color='blue')
                    plt.plot([-Lb, -(l+L)],[maxv2, maxv2], color='blue')
                    #plt.plot(Fyp3[0:I+1])
                    plt.xlabel('Arm location')
                    plt.ylabel('Moment')
                    plt.grid(True)
                    plt.show()
                elif opt==8:
                    j=floor(linspace(1,I+1,10, endpoint=True))
                    for i in j:
                        plt.plot([xcw[i,0], xtip[i,0]],[xcw[i,1], xtip[i,1]], color='blue')
                        plt.plot([xp[i,0], xtip[i,0]],[xp[i,1], xtip[i,1]], color='red')
                        plt.plot([xbarcg[i,0],xbarcg[i,0]],[xbarcg[i,1],xbarcg[i,1]],marker='o',color='green')
                        plt.plot([-l*cos(theta[i]), -l*cos(theta[i])],[0,0],marker='o',color='yellow')

                    plt.grid(True)
                    plt.show()
                elif opt==9:
                    plt.plot(T[0:I+1],Vtot[0:I+1])
                    plt.xlabel('Time (s)')
                    plt.ylabel('Velocity (fps)')
                    plt.grid(True)
                    plt.show()
                elif opt==10:
                    plt.plot(Pos[0:J+1,0],Pos[0:J+1,1])
                    plt.grid(True)
                    plt.show()

                      

    
        else:
            return maxR
    else:
        #print 'Sim Ended'
        return maxR




def func(inputs,x1,s,AL):
    #print '\nx1 = ' + str(x1)

    #print s
    #print x1
    inputs[2]=float(inputs[2]+x1*s[0])
    inputs[3]=float(AL-inputs[2])
    inputs[4]=float(inputs[4]+x1*s[1])
    
    #print 'Trying L = '+ str(inputs[2]) + ' l = ' + str(inputs[3]) + ' d = ' + str(inputs[4])
    if inputs[2]<0.1 or inputs[3]<0.1 or inputs[4]<0:
        out = -(abs(inputs[2])+abs(inputs[3])+abs(inputs[4]))**2
    else:
        out =  main(inputs,runopts)
    #print '*Range = ' + str(out)
    return out



          
def gss(inputs,a,b,eps,N,s,AL):
    
    C=(-1+sqrt(5))/2
    x1=C*a+(1-C)*b
    fx1 = float(func(inputs[:],x1,s,AL))
    x2 = (1-C)*a + C*b
    fx2 = float(func(inputs[:],x2,s,AL))
    
    fa = float(func(inputs[:],a,s,AL))
    fb = float(func(inputs[:],b,s,AL))
    print '---------------------------------------------------------------------------------'
    print '    x1                x2             f(x1)           f(x2)           b - a'
    print '---------------------------------------------------------------------------------'
    print ' %(#)9.5g\t %(#2)9.5g\t %(#3)9.9g\t %(#4)9.9g\t  %(#5)9.9g' %{'#': float(x1),'#2':float(x2),'#3':float(fx1),'#4':float(fx2),'#5':float(b-a)}
    #print str(x1)+'\t'+str(x2)+'\t'+str(fx1)+'\t'+str(fx2)+'\t'+str(x2-x1)
    
    for i in range(0,N):
        #print 'fx1 ' + str(fx1) + ' fx2 ' + str(fx2)
        if fx1 > fx2:
            #print 'fx1 larger'
            fb=fx2
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = C*a + (1-C)*b
            fx1 = float(func(inputs[:],x1,s,AL))
            fc=fx1
            c=x1
        else:
            #print 'fx2 larger'
            fa=fx1
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = (1-C)*a + C*b
            fx2 = float(func(inputs[:],x2,s,AL))
            fc=fx2
            c=x2

            
        top=(c-a)**2.*(fc-fb)-(c-b)**2.*(fc-fa)
        bot=(c-a)*(fc-fb)-(c-b)*(fc-fa)

        d=c-0.5*top/bot
        fd=float(func(inputs[:],d,s,AL))
        
        if d>a and d<x1:
            b=x2
            fb=fx2
            x2=x1
            fx2=fx1
            x1=d
            fx1=fd
        elif d>x1 and d<x2:
            if abs(x1-d)<abs(x2-d):
                b=x2
                x2=d
                fx2=fd
            else:
                a=x1
                x1=d
                fx1=fd
        elif d>x2 and d<b:
            a=x1
            fa=fx1
            x1=x2
            fx1=fx2
            x2=d
            fx2=fd
            
        print ' %(#)9.5g\t %(#2)9.5g\t %(#3)9.9g\t %(#4)9.9g\t  %(#5)9.9g' %{'#': float(x1),'#2':float(x2),'#3':float(fx1),'#4':float(fx2),'#5':float(b-a)}
        if (abs(b-a) < eps):
            print('succeeded after '+str(i)+ ' steps');
            out=(b+a)/2;
            return out


    print('failed requirements after max steps');
        
    


                


#cProfile.run('main()',sort=1, filename="treb_adapt.cprof")

#stats = pstats.Stats("treb_adapt.cprof")
##stats.print_stats()
#stats.strip_dirs().sort_stats('time').print_stats(20)

if __name__ == "__main__":
    start = time.time()
    runopts = opts(0,0,0)
    inputs = getinputfromfile(runopts)

    [M,m,L,l,d,th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]=inputs
   
    x=zeros([1,2])
    x[0,0]=L
    x[0,1]=d
    s=mat([1,2])

    grad=mat([[0.0],[0.0]])
    grad2=mat([[0.0],[0.0]])
    step=0.01
    breaktol=0.5
    B=mat(eye(2))/10.
    AL=3.0
    print '\n\nStarting Optimization'
    for i in range(0,20):
        print '\n\nIteration = ' + str(i+1)
        
        if i==0:
            f1=main(inputs,runopts)
        else:
            f1=f3

        print 'Starting Range = ' + str(f1)
        
        print 'Calculating inital gradient...'
        input0=[M,m,x[0,0]+step,AL-(x[0,0]+step),x[0,1],th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]
        #print input0
        #print step

        grad[0,0]=(main(input0,runopts)-f1)/step
        
        input0=[M,m,x[0,0],AL-x[0,0],x[0,1]+step,th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]
        grad[1,0]=(main(input0,runopts)-f1)/step

        print 'Initial gradient = ' + str(grad)

        s=-B*grad
        print 'GSS'

        a=gss(inputs,-1,1,10**-4.,50,s,AL)

        xold=x
        
        x=x+a*s.T

        print '\nOriginal x' + str(xold)
        print 'Original range' + str(f1)
        print '****New X = ' + str(x)

        
        input0=[M,m,x[0,0],AL-x[0,0],x[0,1],th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]
        #print input0
        print 'Calculating range with new x'
        f3=main(input0,runopts)
        print '****New range = ' + str(f3)
    
        print 'Calculating second gradient'
        input0=[M,m,x[0,0]+step,AL-(x[0,0]+step),x[0,1],th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]
        grad2[0,0]=(main(input0,runopts)-f3)/step
        input0=[M,m,x[0,0],AL-x[0,0],x[0,1]+step,th,mb,Lb,hmax,hmin,Tol,N,thetai,phii,dia,mode]

        grad2[1,0]=(main(input0,runopts)-f3)/step

        q=grad2-grad
        p=x.T-xold.T
       
        
        B=B+float((1+q.T*B*q)/(q.T*p))*(p*p.T)/(p.T*q)-(p*q.T*B+B*q*p.T)/(q.T*p)


        normg2=sqrt(grad2[0,0]**2+grad2[1,0]**2)
        print '*****Norm of second gradient (stopping crit) = ' + str(normg2)

        inputs=input0
        
        i = i + 1

        time.sleep(5)
        if f3-f1<breaktol:
            runopts.verbose=1
            main(inputs,runopts)
            print 'Break'
            print x
            print f3
            elapsed = time.time()-start
            print 'Program run time = %(#)5.3f min' %{'#': float(elapsed/60)}
            sys.exit()

            
        

 
