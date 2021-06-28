      PROGRAM MAIN
 
      implicit none

      double precision Evib,Erot,Norm,k_Bolz,eV2Joule,Tn 
      double precision vel0,vel,dvel,alpha
      double precision RTn_up,RTn_down,Ein
      double precision R_EinTn,RTn,EnergyVel,H2mass

      integer v,j,i,m,n,vib,rot,Tempstep
      integer vmax,jmax
      integer nvelmax

      parameter(vmax=1,jmax=11,Tempstep=8)
      parameter(k_Bolz= 1.38066E-23)  
      parameter(eV2Joule=1.602E-19)
      parameter(H2mass=1.66056E-27*2.0d0)
      parameter(nvelmax=8800,dvel=1.0d0)
      dimension Evib(5),Erot(20),Tn(Tempstep)
      dimension vel0(Tempstep),alpha(Tempstep),EnergyVel(Tempstep)
      dimension RTn_up(Tempstep),RTn_down(Tempstep),RTn(Tempstep)      

C     read data.......
      call  Lecture
C++++++++++++++++++++++++++++

C unit m/s
      vel0(1)=2662.0d0
      vel0(2)=3388.0d0
      vel0(3)=3956.0d0
      vel0(4)=4404.0d0
      vel0(5)=4741.0d0
      vel0(6)=5132.0d0
      vel0(7)=5436.0d0
      vel0(8)=5737.0d0

C alpha unit m/s
      alpha(1)=314.10d0
      alpha(2)=572.20d0
      alpha(3)=795.50d0
      alpha(4)=1007.0d0
      alpha(5)=1165.0d0
      alpha(6)=1346.0d0
      alpha(7)=1528.0d0
      alpha(8)=1705.0d0

C Energy (eV)
      EnergyVel(1)=0.0730d0
      EnergyVel(2)=0.120d0
      EnergyVel(3)=0.160d0 
      EnergyVel(4)=0.210d0
      EnergyVel(5)=0.240d0 
      EnergyVel(6)=0.290d0 
      EnergyVel(7)=0.330d0 
      EnergyVel(8)=0.350d0 


      DO i=1,Tempstep

         Tn(i)=(i-1)*200.0d0+300.0d0
         RTn_up(i)=0.d0
         RTn_down(i)=0.d0

         DO n=1,nvelmax

C           incident energy Ein in eV           
C            print*,"nvelmax", n
            Ein=(1/2.0d0)*H2mass*(dvel*n)**2/eV2Joule
C            print*,"lmax"
            CALL Get_R_EinTn(Ein, Tn(i), R_EinTn)

C      print*,"nvelmax, Ein, Tn(i),R_EinTn",n, Ein, Tn(i),R_EinTn
            RTn_up(i)=RTn_up(i)+((dvel*n)**3 ) *
     &           exp(-(dvel*n-vel0(i))**2/alpha(i)**2)*
     &           R_EinTn*
     &           dvel 

            RTn_down(i)=RTn_down(i)+((dvel*n)**3 ) *
     &           exp(-(dvel*n-vel0(i))**2/alpha(i)**2)*
     &           dvel         
 
         ENDDO
         RTn(i)=RTn_up(i)/RTn_down(i)
      write(*,*) "Temperature(K),Ekin(eV),ReactionProb",
     &                                   Tn(i),EnergyVel(i),RTn(i)
      ENDDO
     

      END   

C===================================================

      SUBROUTINE Get_R_EinTn(Ein, Tn, R_EinTn )

      double precision Evib,Erot,Norm,k_Bolz,eV2Joule,Tn
      double precision vel0,vel,dvel,alpha,velmax,Fb
      double precision R_EinTn,Pr,Ein,ratio_Para_ortho,Ev1j

      integer v,j,i,m,n,vib,rot,Tempstep
      integer vmax,jmax,v1jmax
      integer nvelmax,ParaORortho

      parameter(vmax=2,jmax=11,v1jmax=8)
      parameter(k_Bolz= 1.38066E-23)
      parameter(eV2Joule=1.602E-19)
      parameter(dvel=1.0d0)
      dimension Evib(6), Ev1j(20), Erot(20),Fb(vmax+2,jmax+1)

   
C  Energy unit in eV, vib1=v0 ,rot1=j0
      Evib(1)=0.262770
      Evib(2)=0.766940
    

      Erot(1)=0.262770
      Erot(2)=0.277099 
      Erot(3)=0.305628
      Erot(4)=0.348097
      Erot(5)=0.404130
      Erot(6)=0.473242
      Erot(7)=0.554856
      Erot(8)=0.648312
      Erot(9)=0.752892
      Erot(10)=0.867829
      Erot(11)=0.992326
      Erot(12)=1.125570
      Erot(13)=1.266747
      Erot(14)=1.415048
      Erot(15)=1.569689

     
      Ev1j(1)=0.766940    
      Ev1j(2)=0.780550
      Ev1j(3)=0.807646      
      Ev1j(4)=0.847983
      Ev1j(5)=0.901203
      Ev1j(6)=0.966847
      Ev1j(7)=1.044367
      Ev1j(8)=1.133143 

 
C# part one normalization factor Norm

         Norm=0.d0
      
         v=1
         DO j=1,jmax

C consider para or orth hydrogen in Normal  
            ParaORortho= mod(j,2)
            if( ParaORortho.eq.0 ) then
              ratio_Para_ortho=3.d0
            else
              ratio_Para_ortho=1.d0
            endif
C            print*,"j,ratio :",j-1,ratio_Para_ortho 

            Norm=Norm+exp(-Evib(v)*eV2Joule/(k_Bolz*Tn))* 
     &            (2*j-1)*exp(-(Erot(j)-Erot(1))*
     &            eV2Joule/(0.80d0*k_Bolz*Tn) )*
     &            ratio_Para_ortho    
      
            
         ENDDO
         v=2
         
         DO j=1, v1jmax

           ParaORortho= mod(j,2)
           if( ParaORortho.eq.0 ) then
               ratio_Para_ortho=3.d0
           else
               ratio_Para_ortho=1.d0
           endif
C         print*,"j,ratio :",j-1,ratio_Para_ortho  
           Norm=Norm+exp(-Evib(v)*eV2Joule/(k_Bolz*Tn))*
     &           (2*j-1)*exp(-(Ev1j(j)-Ev1j(1))*
     &           eV2Joule/(0.80d0*k_Bolz*Tn) )*
     &           ratio_Para_ortho
         ENDDO  

C=++++++++++++++++++++++++++++++++++++++++++++++++++
C# part 2, Boltzmann factor Fb
         v=1
         DO j=1,jmax
            ParaORortho= mod(j,2)
            if( ParaORortho.eq.0 ) then
                ratio_Para_ortho=3.d0
            else
                ratio_Para_ortho=1.d0
            endif

            Fb(v,j)=exp(-Evib(v)*eV2Joule/(k_Bolz*Tn))*
     &        ratio_Para_ortho*
     &        (2*j-1)*exp(-(Erot(j)-Erot(1))*
     &        eV2Joule/(0.80d0*k_Bolz*Tn))/Norm
         ENDDO

         v=2
         DO j=1,v1jmax
           ParaORortho= mod(j,2)
           if( ParaORortho.eq.0 ) then
               ratio_Para_ortho=3.d0
           else
               ratio_Para_ortho=1.d0
           endif
           Fb(v,j)=exp(-Evib(v)*eV2Joule/(k_Bolz*Tn))*
     &        ratio_Para_ortho*
     &        (2*j-1)*exp(-(Ev1j(j)-Ev1j(1) )*
     &        eV2Joule/(0.80d0*k_Bolz*Tn))/Norm
        ENDDO

C        print*,"v1jmax= ",v1jmax
C+++++++++++++++++++++++++++++++++++++++++++++++++++++
C# part 3, energy resolved reaction prob: R(Ein; Tn)
         v=1
         R_EinTn=0.d0
         DO j=1,jmax
             call GetQCT_reactionProb(v,j,Ein,Pr)      
             R_EinTn=R_EinTn+Fb(v,j)*Pr
         ENDDO

         v=2
         j=1
         call GetQCT_reactionProb(v,j,Ein,Pr)
         R_EinTn=R_EinTn+Fb(v,j)*Pr


         END




      SUBROUTINE Lecture
      
      implicit none
      integer j,mj,nE,k,m,n,Jmax,mjmax,nEmax
      double precision P,a,b,c,averP  
 
      parameter(Jmax=12,mjmax=12,nEmax=60)
      dimension P(Jmax+12,Jmax+1,nEmax),averP(2,mjmax,nEmax)
 
      COMMON/ProbDATA/averP

      OPEN(10,FILE='j0mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(0+1,0+1,n)
      ENDDO
      CLOSE(10)
      
      OPEN(10,FILE='j1mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(1+1,0+1,n)
      ENDDO
      CLOSE(10) 

      OPEN(10,FILE='j1mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(1+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j2mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(2+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j2mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(2+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j2mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(2+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j3mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(3+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j3mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(3+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j3mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(3+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j3mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(3+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j4mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(4+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j4mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(4+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j4mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(4+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j4mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(4+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j4mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(4+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j5mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(5+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j6mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(6+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j7mj7_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(7+1,7+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj7_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,7+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j8mj8_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(8+1,8+1,n)
      ENDDO
      CLOSE(10)

C    j=9
      OPEN(10,FILE='j9mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj7_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,7+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj8_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,8+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j9mj9_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(9+1,9+1,n)
      ENDDO
      CLOSE(10)

C    j=10

      OPEN(10,FILE='j10mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,0+1,n)
      ENDDO
      CLOSE(10)


      OPEN(10,FILE='j10mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj7_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,7+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj8_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,8+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj9_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,9+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='j10mj10_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(10+1,10+1,n)
      ENDDO
      CLOSE(10)

C++++++++++++++++++++++++ v=1 states j=0-7 +++++++++++++++++++++++++++
C  v1j0mj0
      OPEN(10,FILE='v1j0mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(11+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j1mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(12+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j1mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(12+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j2mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(13+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j2mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(13+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j2mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(13+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j3mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(14+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j3mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(14+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j3mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(14+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j3mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(14+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j4mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(15+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j4mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(15+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j4mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(15+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j4mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(15+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j4mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(15+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j5mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(16+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,0+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,1+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,2+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,3+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,4+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,5+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j6mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(17+1,6+1,n)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='v1j7mj0_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,0+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj1_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,1+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj2_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,2+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj3_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,3+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj4_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,4+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj5_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,5+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj6_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,6+1,n)
      ENDDO
      CLOSE(10)
      OPEN(10,FILE='v1j7mj7_prob.out',STATUS='OLD')
      DO n=1,nEmax
        READ(10,*) a,b,c,P(18+1,7+1,n)
      ENDDO
      CLOSE(10)










C Average prob
C  use cycle

c---old average 
      DO n=1,nEmax
        averP(1,1,n)=P(1,1,n)
          write(50,*) "j0 average",n*0.027211d0,averP(1,1,n)

        averP(1,2,n)=(P(2,1,n)+P(2,2,n)*2)/3.d0
          write(51,*) "j1 average",n*0.027211d0,averP(1,2,n)

        averP(1,3,n)=(P(3,1,n)+P(3,2,n)*2+P(3,3,n)*2)/5.d0
          write(52,*) "j2 average",n*0.027211d0,averP(1,3,n)

        averP(1,4,n)=(P(4,1,n)+P(4,2,n)*2+P(4,3,n)*2+P(4,4,n)*2)/7.d0
          write(53,*) "j3 average",n*0.027211d0,averP(1,4,n)

        averP(1,5,n)=(P(5,1,n)+P(5,2,n)*2+P(5,3,n)*2+P(5,4,n)*2+
     &              P(5,5,n)*2)/9.d0
          write(54,*) "j4 average",n*0.027211d0,averP(1,5,n)  

        averP(1,6,n)=(P(6,1,n)+P(6,2,n)*2+P(6,3,n)*2+P(6,4,n)*2+
     &              P(6,5,n)*2+ P(6,6,n)*2)/11.d0
          write(55,*) "j5 average",n*0.027211d0,averP(1,6,n)  

        averP(1,7,n)=(P(7,1,n)+P(7,2,n)*2+P(7,3,n)*2+P(7,4,n)*2+
     &              P(7,5,n)*2+ P(7,6,n)*2+P(7,7,n)*2 )/13.d0
          write(56,*) "j6 average",n*0.027211d0,averP(1,7,n)

        averP(1,8,n)=(P(8,1,n)+P(8,2,n)*2+P(8,3,n)*2+P(8,4,n)*2+
     &              P(8,5,n)*2+ P(8,6,n)*2+P(8,7,n)*2 +P(8,8,n)*2)/15.d0
          write(57,*) "j7 average",n*0.027211d0,averP(1,8,n)      

        averP(1,9,n)=(P(9,1,n)+P(9,2,n)*2+P(9,3,n)*2+P(9,4,n)*2+
     &              P(9,5,n)*2+ P(9,6,n)*2+P(9,7,n)*2 +P(9,8,n)*2+
     &              P(9,9,n)*2 )/17.d0
          write(58,*) "j8 average",n*0.027211d0,averP(1,9,n)

        averP(1,10,n)=(P(10,1,n)+P(10,2,n)*2+P(10,3,n)*2+P(10,4,n)*2+
     &             P(10,5,n)*2+ P(10,6,n)*2+P(10,7,n)*2 +P(10,8,n)*2+
     &             P(10,9,n)*2+ P(10,10,n)*2 )/19.d0
          write(59,*) "j9 average",n*0.027211d0,averP(1,10,n)

        averP(1,11,n)=(P(11,1,n)+P(11,2,n)*2+P(11,3,n)*2+P(11,4,n)*2+
     &             P(11,5,n)*2+ P(11,6,n)*2+P(11,7,n)*2+P(11,8,n)*2+
     &             P(11,9,n)*2+ P(11,10,n)*2 +P(11,11,n)*2 )/21.d0
          write(60,*) "j10 average",n*0.027211d0,averP(1,11,n)


  
        averP(2,1,n)=P(12,1,n)
          write(70,*) "v1j0 average",n*0.027211d0,averP(2,1,n)
    
        averp(2,2,n)=(P(13,1,n)+P(13,2,n)*2 )/3.0d0
          write(71,*) "v1j1 average",n*0.027211d0,averP(2,2,n)

        averp(2,3,n)=(P(14,1,n)+P(14,2,n)*2 +P(14,3,n)*2)/5.0d0
          write(72,*) "v1j2 average",n*0.027211d0,averP(2,3,n)
      
        averp(2,4,n)=(P(15,1,n)+P(15,2,n)*2 +P(15,3,n)*2 +P(15,4,n)*2)/7.0d0
          write(73,*) "v1j3 average",n*0.027211d0,averP(2,4,n)

        averp(2,5,n)=(P(16,1,n)+P(16,2,n)*2 +P(16,3,n)*2 +P(16,4,n)*2
     &               +P(16,5,n)*2 )/9.0d0
          write(74,*) "v1j4 average",n*0.027211d0,averP(2,5,n)

        averp(2,6,n)=(P(17,1,n)+P(17,2,n)*2 +P(17,3,n)*2 +P(17,4,n)*2
     &               +P(17,5,n)*2 +P(17,6,n)*2 )/11.0d0
          write(75,*) "v1j5 average",n*0.027211d0,averP(2,6,n)

        averp(2,7,n)=(P(18,1,n)+P(18,2,n)*2 +P(18,3,n)*2 +P(18,4,n)*2
     &               +P(18,5,n)*2 +P(18,6,n)*2 + P(18,7,n)*2  )/13.0d0
          write(76,*) "v1j6 average",n*0.027211d0,averP(2,7,n)

        averp(2,8,n)=(P(19,1,n)+P(19,2,n)*2 +P(19,3,n)*2 +P(19,4,n)*2
     &               +P(19,5,n)*2 +P(19,6,n)*2 + P(19,7,n)*2 
     &                +P(19,8,n)*2 )/15.0d0
          write(77,*) "v1j7 average",n*0.027211d0,averP(2,8,n)

      ENDDO

      print*,'finished average'


      END
   

      SUBROUTINE GetQCT_reactionProb(vib,rot,Ein,Prob)
      IMPLICIT NONE

      integer i,n,j,k,l,jmax,vib,rot,nEmax,ParaORortho
      double precision x,y,vv,vvv,Ein,Prob
      double precision b,c,d,s,r,rx,averP
      double precision eV2Joule,ratio_Para_ortho      


 
      parameter(nEmax=60,jmax=12,n=60)
      parameter(eV2Joule=1.602E-19)

      dimension x(n),y(n),vv(n),vvv(n)
      dimension b(n),c(n),d(n),averP(2,12,nEmax)
      
      COMMON/ProbDATA/averP


      DO i=1,nEmax
         x(i)=i*0.027211d0     
         y(i)=averP(vib,rot,i)
      ENDDO
  
      call spline(nEmax, x, y, b, c, d)
C      DO i=1,n
C         write(25,*)  b(i), c(i), d(i)
C      ENDDO        

C     Ein are in eV already
C      Ein=2.0d0*(1.66056E-27)*(dvel*n)**2/eV2Joule
      rx=Ein
C
      k=int(Ein/0.027211d0)
      s=y(k)+b(k)*(rx-x(k))+c(k)*(rx-x(k))**2 +d(i)*(rx-x(k))**3
      Prob=s
C      Prob=s*ratio_Para_ortho

C# Check????????????????????????????????????????????????????????????????
C????????????????    ?????
C      DO j=1, jmax 
C         r=((j-1)*(n-1)/(jmax-1) +1)
C         rx=((j-1)*(n-1)/(jmax-1) +1)*0.027211d0
C         k=int(r)
C         s=y(k)+b(k)*(rx-x(k))+c(k)*(rx-x(k))**2 +d(i)*(rx-x(k))**3
C         write(26,*) rx,s
C      ENDDO
C  ????????????????????????????????
 

      END



      subroutine spline(n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end


