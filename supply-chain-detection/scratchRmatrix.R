library(plot3D)


#1-eps,eps case

eps=0.1


A = matrix(c(1-eps,eps,0,0,
             1-eps,0,eps,0,
             1-eps,0,0,eps,
             eps,1-eps,0,0,
             0,1-eps,eps,0,
             0,1-eps,0,eps,
             eps,0,1-eps,0,
             0,eps,1-eps,0,
             0,0,1-eps,eps,
             eps,0,0,1-eps,
             0,eps,0,1-eps,
             0,0,eps,1-eps),
           nrow=12,ncol=4,byrow=TRUE)


t(A)%*%A
solve(t(A)%*%A)%*%t(A)


n=12
m=4
e=eps
denom = (n/m)*(1-(2*e)+2*(e^2))


epsVec = seq(0.01,0.99,by = 0.01)
mVec = seq(5,30,by = 1)

c=10
zPrim = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=m*(m-1)*c
    e = epsVec[j]
    zPrim[i,j]=((1-e)*(m-1)-2*e+2*(e^2))/(n*(((m-1)/m)-2*e+2*e^2))
  }
}
contour(mVec,epsVec,zPrim,main='Primary children entries WRT m, epsilon',xlab='m',ylab='epsilon') 
persp(mVec,epsVec,zPrim,main='Primary children entries WRT m, epsilon',xlab='m',ylab='epsilon',d=1,theta=100,phi=30,ticktype='detailed') 

zSecond = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=m*(m-1)*c
    e = epsVec[j]
    zSecond[i,j]=((e)*(m-1)-2*e+2*(e^2))/(n*(((m-1)/m)-2*e+2*e^2))
  }
}
contour(mVec,epsVec,zSecond,main='Secondary children entries WRT m, epsilon',xlab='m',ylab='epsilon') 
persp(mVec,epsVec,zSecond,main='Secondary children entries WRT m, epsilon',xlab='m',ylab='epsilon',d=1,theta=100,phi=30,ticktype='detailed') 

zNone = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=m*(m-1)*c
    e = epsVec[j]
    zNone[i,j]=(-2*e+2*(e^2))/(n*(((m-1)/m)-2*e+2*e^2))
  }
}
contour(mVec,epsVec,zNone,main='Non-children entries WRT m, epsilon',xlab='m',ylab='epsilon') 
persp(mVec,epsVec,zNone,main='Non-children entries WRT m, epsilon',xlab='m',ylab='epsilon',d=1,theta=100,phi=30,ticktype='detailed') 

#######################################################
#1-eps,eps/(m-1) case
eps=0.1
m=4
e_m = eps/(m-1)

A = matrix(c(1-eps,e_m,e_m,e_m,
             1-eps,e_m,e_m,e_m,
             1-eps,e_m,e_m,e_m,
             e_m,1-eps,e_m,e_m,
             e_m,1-eps,e_m,e_m,
             e_m,1-eps,e_m,e_m,
             e_m,e_m,1-eps,e_m,
             e_m,e_m,1-eps,e_m,
             e_m,e_m,1-eps,e_m,
             e_m,e_m,e_m,1-eps,
             e_m,e_m,e_m,1-eps,
             e_m,e_m,e_m,1-eps),
           nrow=12,ncol=4,byrow=TRUE)

t(A)%*%A
solve(t(A)%*%A)%*%t(A)

n=12
m=4
(n/m)*((1-eps)^2 + (eps^2)/(m-1))
(n/m)*(2*eps*(1-eps)/(m-1) + (m-2)*(eps^2)/((m-1)^2))

#put contours
c=10
#diagonals of (t(A)*A)^-1
  #detailed equation
zPrim = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=c*m
    e = epsVec[j]
    if ((e >= (m-1)/m-0.05) & (e <= (m-1)/m+0.05)) { zPrim[i,j] = NA }
    else {
      numer = 1 - (2/(m-1))*e + m*(e^2)/((m-1)^2)
      denom = (n/m)*(1 - 2*m*e/(m-1) + (m^2)*(e^2)/((m-1)^2)) 
      zPrim[i,j]= numer/denom
    }
  }
}
persp3D(mVec,epsVec,zPrim,theta=-60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="epsilon", zlab="diagonal value", 
        main="AtA^(-1) diagonals")

#non-diagonals of (t(A)*A)^-1
zPrim = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=c*m
    e = epsVec[j]
    if ((e >= (m-1)/m-0.05) & (e <= (m-1)/m+0.05)) { zPrim[i,j] = NA }
    else {
      numer =  (-2/(m-1))*e + m*(e^2)/((m-1)^2)
      denom = (n/m)*(1 - 2*m*e/(m-1) + (m^2)*(e^2)/((m-1)^2)) 
      zPrim[i,j]= numer/denom
    }
  }
}
persp3D(mVec,epsVec,zPrim,theta=-60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="epsilon", zlab="non-diagonal value", 
        main="AtA^(-1) non-diagonals")

#primary children of (t(A)*A)^-1*t(A)
zPrim = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=c*m
    e = epsVec[j]
    if ((e >= (m-1)/m-0.05) & (e <= (m-1)/m+0.05)) { zPrim[i,j] = NA }
    else {
      numer =  1-e - (2/(m-1))*e + m*(e^2)/((m-1)^2)
      denom = (n/m)*(1 - 2*m*e/(m-1) + (m^2)*(e^2)/((m-1)^2)) 
      zPrim[i,j]= numer/denom
    }
  }
}
persp3D(mVec,epsVec,zPrim,theta=-60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="epsilon", zlab="entry value", 
        main="(AtA^(-1))*At primary children")

#secondary children of (t(A)*A)^-1*t(A)
zPrim = matrix(0,nrow=length(mVec),ncol=length(epsVec))
for (i in 1:length(mVec)){
  for (j in 1:length(epsVec)){
    m = mVec[i]
    n=c*m
    e = epsVec[j]
    if ((e >= (m-1)/m-0.05) & (e <= (m-1)/m+0.05)) { zPrim[i,j] = NA }
    else {
      numer =  (-1/(m-1))*e + m*(e^2)/((m-1)^2)
      denom = (n/m)*(1 - 2*m*e/(m-1) + (m^2)*(e^2)/((m-1)^2)) 
      zPrim[i,j]= numer/denom
    }
  }
}
persp3D(mVec,epsVec,zPrim,theta=-30, phi=-5, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="epsilon", zlab="entry value", 
        main="(AtA^(-1))*At secondary children")




#######################################################
#rho-decay case
rho = 0.3
m=4
n=m*(m-1)
c = ((1-rho)/(1-(rho^m)))
rVec = c(rho^0,rho^1,rho^2,rho^3)*c
A = matrix(c(rVec[1],rVec[2],rVec[3],rVec[4],
             rVec[1],rVec[3],rVec[4],rVec[2],
             rVec[1],rVec[4],rVec[2],rVec[3],
             rVec[2],rVec[1],rVec[4],rVec[3],
             rVec[3],rVec[1],rVec[2],rVec[4],
             rVec[4],rVec[1],rVec[3],rVec[2],
             rVec[3],rVec[4],rVec[1],rVec[2],
             rVec[4],rVec[2],rVec[1],rVec[3],
             rVec[2],rVec[3],rVec[1],rVec[4],
             rVec[4],rVec[3],rVec[2],rVec[1],
             rVec[2],rVec[4],rVec[3],rVec[1],
             rVec[3],rVec[2],rVec[4],rVec[1]),
           nrow=12,ncol=4,byrow=TRUE)



ata = t(A)%*%A
ata_inv = solve(ata)
proj_mat = ata_inv%*%t(A)

### CHECKING EQUATIONS
#ata calculations
S1 = 0
S2 = 0
for (j in 0:(m-1)){
  S1 = S1 + rho^(2*j)
}
for (i in 0:(m-1)){
  for (j in 0:(m-1)){
    if (j>i){
    
      S2 = S2 + rho^(i+j)  
    }
  }
}
S2 = 2*S2

ata[1,1]
(c^2)*(n/m)*S1
ata[1,2]
(c^2)*S2
### ata WORKS

### ata_inv calcs
a = (c^2)*(n/m)*S1
b = (c^2)*S2

Dnum = (m-2)*S2/(m-1) + S1
NDnum = -1*S2/(m-1)
denom = ((m-1)*S1 - S2)*(S1+S2)
c_1 = c^(-2)

ata_inv[1,1]
inv_D = Dnum*c_1/denom
ata_inv[1,2]
inv_ND = NDnum*c_1/denom

### ata_inv WORKS

### proj_mat calcs
l=0
numer = ((rho^l)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
numer/denom
proj_mat
### proj_mat WORKS


# S1, S2 WRT m, rho
mVec = seq(5,30,by=1)
rhoVec = seq(0.01,0.99,by=0.01)

zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    zPrim[i,j]= S1
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=15, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="S_1", 
        main="S1 values WRT m, rho")

zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    zPrim[i,j]= 2*S2
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="S_2", 
        main="S_2 values WRT m, rho")

#denominator WRT m,rho
zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    denom = ((m-1)*S1 - S2)*(S1+S2)
    zPrim[i,j]= denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=30,phi=25, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="Denominator value", 
        main="Denominator term WRT m, rho")

plot(mVec,zPrim[,99]/mVec,xlab='m',ylab='S2/(S1*m)',type='l',lwd=4,main='S2/(S1*m) for rho=0.99',col='orangered')

# (AtA)^(-1) WRT m, rho
mVec = seq(5,30,by=1)
rhoVec = seq(0.01,0.9,by=0.01)

zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    Dnum = (m-2)*S2/(m-1) + S1
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    
    zPrim[i,j] = Dnum*c_1/denom
    }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=15, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="diagonal", 
        main="(AtA)^(-1) diagonals WRT m, rho")


zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    NDnum = -1*S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = NDnum*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="non-diagonal", 
        main="(AtA)^(-1) non-diagonals WRT m, rho")





# (AtA)^(-1)*At entries for l=0,1,2,5,m
mVec = seq(4,30,by=1)
rhoVec = seq(0.01,0.9,by=0.01)

l = 0
zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    numer = ((rho^l)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
    
    numer = (rho^l)*c*(S1+S2) - S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = numer*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=15, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 0")

persp3D(mVec[5:27],rhoVec[1:70],zPrim[5:27,1:70],theta=160, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 0")


l = 1
zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    numer = ((rho^l)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
    
    numer = (rho^l)*c*(S1+S2) - S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = numer*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 1")

persp3D(mVec[10:27],rhoVec[1:80],zPrim[10:27,1:80],theta=160, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 1")

l = 2
zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    numer = ((rho^l)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
    
    numer = (rho^l)*c*(S1+S2) - S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = numer*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 2")

persp3D(mVec[10:27],rhoVec[1:80],zPrim[10:27,1:80],theta=160, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 2")

l = 3
zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    numer = ((rho^l)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
    
    numer = (rho^l)*c*(S1+S2) - S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = numer*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 3")

persp3D(mVec[10:27],rhoVec[1:80],zPrim[10:27,1:80],theta=160, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = 3")


zPrim = matrix(0,nrow=length(mVec),ncol=length(rhoVec))
for (i in 1:length(mVec)){
  for (j in 1:length(rhoVec)){
    m = mVec[i]
    r = rhoVec[j]
    n = m*(m-1)
    c = ((1-r)/(1-(r^m))) 
    
    S1 = 0
    for (k in 0:(m-1)){
      S1 = S1 + r^(2*k)
    }
    S2 = 0
    for (ii in 0:(m-1)){
      for (jj in 0:(m-1)){
        if (jj>ii){
          S2 = S2 + r^(ii+jj)  
        }
      }
    }
    S2=2*S2
    
    numer = ((rho^m)*(c^(-1))*(S1+S2) - (c^(-2))*S2*(1/(m-1)))
    
    numer = (rho^l)*c*(S1+S2) - S2/(m-1)
    denom = ((m-1)*S1 - S2)*(S1+S2)
    c_1 = c^(-2)
    zPrim[i,j] = numer*c_1/denom
  }
}
persp3D(mVec,rhoVec,zPrim,theta=60, phi=10, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = m")

persp3D(mVec[10:27],rhoVec[1:80],zPrim[10:27,1:80],theta=160, phi=20, axes=TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed", xlab="m", ylab="rho", zlab="entry value", 
        main="(AtA)^(-1)*At entries WRT m, rho: l = m")




















#######################################################



A = matrix(c(1-eps,1-eps,1-eps,eps,0,eps,0,eps,eps,0,
             eps,eps,0,1-eps,1-eps,1-eps,1-eps,0,0,eps,
             0,0,eps,0,eps,0,eps,1-eps,1-eps,1-eps),
           nrow=10,ncol=3)

AtA = t(A)%*%A

AtA_inv = solve(AtA)

alpha=diag(3)

#1
newU = c(1-eps,eps,0)
v = newU %o% newU
alpha=alpha+v
#2
newU = c(1-eps,eps,0)
v = newU %o% newU
alpha=alpha+v
#3
newU = c(1-eps,0,eps)
v = newU %o% newU
alpha=alpha+v
#4
newU = c(eps,1-eps,0)
v = newU %o% newU
alpha=alpha+v
#5
newU = c(0,1-eps,eps)
v = newU %o% newU
alpha=alpha+v
#6
newU = c(eps,1-eps,0)
v = newU %o% newU
alpha=alpha+v
#7
newU = c(0,1-eps,eps)
v = newU %o% newU
alpha=alpha+v
#8
newU = c(eps,0,1-eps)
v = newU %o% newU
alpha=alpha+v
#9
newU = c(eps,0,1-eps)
v = newU %o% newU
alpha=alpha+v
#10
newU = c(0,eps,1-eps)
v = newU %o% newU
alpha=alpha+v

alpha = alpha-diag(3)











