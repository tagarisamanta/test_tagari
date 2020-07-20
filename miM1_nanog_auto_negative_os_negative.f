c       GILLESPIE SIMULATION
C     DATE: 8th AUGUST, 2014
c.............................................................................
c     USING "DIRECT METHOD" OF  
c     GILLESPIE'S "STOCHASTIC SIMULATION ALGORITHM".
c     dated 8.08.2014
c     Author: Tagari Samanta
c
      implicit real*8 (a-h,q-z)
       Real pp1

      parameter(nu=54)
      parameter(nn=2000)
      dimension a_nu(nu),x11(nn),y11(nn),x12(nn),y12(nn),sqy1(nn)
      dimension multi1(nn),diff1(nn),sq1(nn),sqy2(nn),sqx2(nn),sqt2(nn)
      dimension multi2(nn),diff2(nn),sq2(nn),sqt1(nn),sqx1(nn)
      dimension x13(nn), y13(nn)
            
c     a_nu = h_nu*c_nu
c     h_nu = number of distinct combinations of reactant molecules for 
c          reaction R_nu that are present in V at time t.
c     c_nu = reaction parameter (k_nu = V*c_nu)
c.......................................................................
c      open(50,file='all.out',status='unknown')
c      open(2401,file='Ga.out',status='unknown')
c      open(2402,file='Ga1.out',status='unknown')
c      open(2403,file='Ga2.out',status='unknown')
c      open(2404,file='protein_P1.out',status='unknown') 
c      open(2405,file='protein_P2.out',status='unknown')
c     open(2406,file='Total_mRNA_M.out',status='unknown')
c      open(2407,file='oct4_sox2_complex_OS.out',status='unknown')
c      open(2408,file='Oct4.out',status='unknown')
c      open(2409,file='Sox2.out',status='unknown')
c      open(2410,file='mRNA_Oct4.out',status='unknown')
c      open(2411,file='mRNA_Sox2.out',status='unknown')
c      open(2412,file='Total_protein_X.out',status='unknown')
c      open(2413,file='mRNA_M1.out',status='unknown')
c      open(2414,file='mRNA_M2.out',status='unknown')
c      open(2415,file='micro_RNA_miR.out',status='unknown')
c      open(2416,file='microRNA_M1_complex_miM1.out',status='unknown')
c      open(2417,file='microRNA_M2_complex_miM2.out',status='unknown')
      open(2501,file='Protein_P1_stat.out',status='unknown') 
      open(2502,file='Protein_P2_stat.out',status='unknown')
      open(2503,file='mRNA1_M1_stat.out',status='unknown')             
      open(2504,file='oct4_sox2_complex_OS_out.out',status='unknown')
      open(2505,file='Oct4_stat.out',status='unknown')
      open(2506,file='Sox2_stat.out',status='unknown')
      open(2507,file='Total_protein_X_stat.out',status='unknown')
      open(2508,file='mRNA2_M2_stat.out',status='unknown')
      open(2509,file='mRNA_M1_M2_stat.out',status='unknown')
      open(2510,file='mRNA_M1_M2_divided.out',status='unknown')
      open(2511,file='mRNA_total.out',status='unknown')
      open(2512,file='m1_m2_stat.out',status='unknown')
      open(2513,file='microRNA_stat.out',status='unknown')
      open(2514,file='microRNA_M1_stat.out',status='unknown')
      open(2515,file='microRNA_M2_stat.out',status='unknown')
      open(2516,file='free_mRNA_and_complex.out',status='unknown')
      open(2517,file='protein_and_dimer.out',status='unknown')
      open(5001,file='mRNA_noise_data.out',status='unknown')
      open(5002,file='protein_noise_data.out',status='unknown')
c      open(5003,file='data.out',status='unknown')
c      open(3001,file='asyn_stat.out',status='unknown')
c      open(3002,file='rest_time.out',status='unknown')
c      open(3003,file='cell11.out',status='unknown')
c      open(3004,file='cell12.out',status='unknown')
c      open(3005,file='cell21.out',status='unknown')
c      open(3006,file='cell22.out',status='unknown')
c      open(3007,file='cell13.out',status='unknown')
c      open(3008,file='cell23.out',status='unknown')
c      open(1000,file='points_noise.out',status='unknown')
c      open(2100,file='phase.out',status='unknown')
c      open(2200,file='G1_phase.out',status='unknown')
c      open(2300,file='S_phase.out',status='unknown')
c      open(2400,file='G2_phase.out',status='unknown')
c      open(2500,file='M_phase.out',status='unknown')
c      open(3000,file='phase_distribution.out',status='unknown')
c      open(10,file='gene_activation_rate.out',status='unknown') 
c      open(12,file='transcription_rate_values.out',status='unknown') 
c    ------------------------------------------------------------------
c     STEP 0: INITIALIZATION
c     -----------------------------------------------------------------
        nk=2000
         aobs_time=240.0
         ccount=0.0
         bcount=0.0
          sigma3=192.0
c        idum=0
        l100=1
c.....................................................................
c      Parameter values for the problem
c.....................................................................
       ak1=0.663
       ak2=0.609
       ak3=1.0
       ak4=0.00013
       ak5=50.0
       ak6=0.6
       ak7=0.0034
       ak8=0.001
       ak9=0.0001
       ak10=0.46
       ak11=0.005
       ak12=0.0
       ak13=0.0001
       ak14=1.0
       ak15=0.00195
       ak16=0.0104
       ak17=0.0104
       ak18=0.01
       ak19=0.01
       ak20=0.1
       ak21=60.0
       ak22=0.00495
       ak23=0.009
       ak24=0.02
       ak25=0.0005
       ak26=1.0
       ak27=0.00023
       ak28=0.01
       ak29=0.0005
       ak30=1.0
       ak31=0.11
       ak32=0.0022
       ak33=0.66
       ak34=0.01
       ak35=0.055
       ak36=0.00277
       ak37=0.05
       ak38=0.01
       ak39=0.0005
       ak40=8.0
       ak81=0.0085
       ak82=0.0125
       ak60=1.0
       ak401=1.0
       ak331=1.0
       ak341=0.78
       ak351=2.0
       ak361=0.55
       ak371=0.3
       ak381=0.0
       ak382=0.0
       ak390=0.2
       ak395=0.2
       ak391=2.3
       ak392=1.7
       ak393=0.4
       ak394=0.103
       ak47=1.0
       ak48=0.3
       ak49=10.0
       ak50=1.46
       ak57=1.0
       ak58=0.35
       ak59=10.0
       ak601=1.46 
       ak611=3.0
       ak1011=3.0
       ak2311=0.00000167
       ak2312=10.0

c...................................................................
       ak90=((ak401+(ak331*ak381**ak351)/(ak341**ak351
     +    +ak381**ak351))*(1-(ak361*ak382)
     +    +(ak371*ak382**2))) 
       ak61=(ak6+(ak391*ak381))*(ak47+((ak48-ak47)*ak382**ak49
     +    /(ak50**ak49+ak382**ak49))) 
       ak101=(ak10+(ak392*ak381))*(ak57+((ak58-ak57)*ak382**ak59
     +     /(ak601**ak59+ak382**ak59)))      
       write(*,*) ak90, ak61, ak101        
c......................................................................
c        write(*,*)'ccount=',ccount,'bcount=',bcount
       do 300 jj=1,nk
       write(*,*) jj
c      call random_seed
        count=0.0 
        l100=l100+1
        iseed=l100
        useed=dustar(iseed)
        api1 = 4.0*atan(1.0)
c......................................................................
c To simulate the 2 cells which have generated due to a division event
c......................................................................
 355         if(ccount.ge.1.0)then
           ti=0.0
            t_f=rest_time+1.0
            end_time=rest_time
              bcount=bcount+1.0
c               write(*,*)'bcount=',bcount,'ti=',ti,'t_f=',t_f
               if(bcount.eq.1.0)then
c......................................................................
c  Read the intial conditions for the 1st daughter cell 
c......................................................................
                   Ga=0.0
                   GaO=0.0
                   GaS=0.0
c                  read(*,3003)a11,a12,a13,a14
                   xP1=xP11
                    zP2=zP21
			  DP1=DP11
			   DP2=DP21
                     yM1=yM51
                     yM2=yM52
                      bOS=bOS1
c               write(*,*)xP,zP2,yM1,yM2,bOS
c                read(*,3004)a15,a16,a17,a18
                     ymO=ymO1
                      ymS=ymS1
                       bOct4=bOct41
                        Sox2=Sox21
c           write(*,*)ymO,ymS,bOct4,Sox2
					ymiR=ymiR1
					 ymiM1=ymiM11
					  ymiM2=ymiM21
c				write(*,*)ymiR,ymiM1,ymiM2					
c.................................................................
c.....................................................................
c Creating cell cycle periods with 16 hr mean and 10% CV for 1st cell
c.....................................................................
                     bb11 = duni()
                        dd11 = duni()
                     api1 = 4.0*atan(1.0)
c
                     if(bb11.le.0.0)then
                   bb11 = 0.00001
                        endif
c
                  ba1  = -2.0*sigma3*sigma3*log(bb11)
                   ba11 = dcos(2.0*api*dd11)
                      sf1  = sqrt(ba1)*ba11
                   cellper = 960.0 + sf1
c                 write(3001,*)cellper,rest_time,'r'
                        else
c.................................................................
c  Read the intial conditions for the 1st daughter cell
c................................................................. 
                   Ga=0.0
                   GaO=0.0
                   GaS=0.0
c                     read(*,3005)b11,b12,b13,b14
                        xP1=xP12
                         zP2=zP22
				 DP1=DP12
				  DP2=DP22
                          yM1=yM61
                          yM2=yM62
                           bOS=bOS2
c                     write(*,*)xP,zP2,yM1,yM2,bOS
c                       read(*,3006)b15,b16,b17,b18
                          ymO=ymO2
                           ymS=ymS2
                            bOct4=bOct42
                             Sox2=Sox22
c                       write(*,*)ymO,ymS,bOct4,Sox2
						ymiR=ymiR2
						 ymiM1=ymiM12
						  ymiM2=ymiM22
c				write(*,*)ymiR,ymiM1,ymiM2					
c                          write(*,*)'ti=',ti,'t_f=',t_f
c.....................................................................
c Creating cell cycle periods with 16 hr mean and 10% CV for 2nd cell
c.....................................................................
                            bb12 = duni()
                        dd12 = duni()
                         api2 = 4.0*atan(1.0)
c
                        if(bb12.le.0.0)then
                             bb12 = 0.00001
                          endif
c
                     ba2  = -2.0*sigma3*sigma3*log(bb12)
                      ba12 = dcos(2.0*api*dd12)
                       sf1  = sqrt(ba2)*ba12
                     cellper = 960.0 + sf1
c                      write(3001,*)cellper,rest_time,'r'
c.................................................................
c..................................................................
c                  Period define                 
c..................................................................
            G1_time=0.25*cellper
            S_time=0.75*cellper
            G2_time=0.85*cellper
        if(diff_time.le.aobs_time) then    
         if(rest_time.le.G1_time) then
             ak111=ak1*2
             ak41=ak4*2
             ak251=ak25*2
             ak291=ak29*2
             ak231=ak23*2
             ak271=ak27*2
             ak62=ak61*1
             ak102=ak101*1
c       write(2100,*) cellper,rest_time,'r', 'G1', ak62, ak102 
c             write(2200,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
            else
           if((rest_time.gt.G1_time) .and. (rest_time.le.S_time)) then
               ak111=ak1*2
               ak41=ak4*2
               ak231=ak23*2
               ak251=ak25*2
               ak271=ak27*2
               ak291=ak29*2
               ak62=ak61*2
               ak102=ak101*2
c       write(2100,*) cellper,rest_time,'r', 'S', ak62, ak102 
c             write(2300,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
              else
            if((rest_time.gt.S_time) .and. (rest_time.le.G2_time)) then
                ak111=ak1
                ak231=ak23
                ak41=ak4
                ak291=ak29
                ak251=ak25
                ak271=ak27
                ak62=ak61*2
                ak102=ak101*2
c       write(2100,*) cellper,rest_time,'r', 'G2', ak62, ak102 
c             write(2400,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
                   else
                    ak111=ak1
                    ak41=ak4
                    ak231=ak23
                    ak251=ak25
                    ak271=ak27
                    ak291=ak29
                    ak62=ak61*2
                    ak102=ak101*2
c      write(2100,*) cellper,rest_time,'r', 'M', ak62, ak102 
c             write(2500,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
               endif     
           endif
        endif 
       endif               

                          endif
c...................................................................... 
c.................................................................                         
c.................................................................
c To simulate the next cell with or without division
c.................................................................
                 else
c.................................................................
c                 write(*,*)'ccount=',ccount,'bcount=',bcount
c.....................................................................
c Creating cell cycle periods with 16 hr mean and 10% CV
c.....................................................................
         bb1 = duni()
            dd1 = duni()
          api = 4.0*atan(1.0)
c
         if(bb1.le.0.0)then
            bb1 = 0.00001
         endif
c
           ba  = -2.0*sigma3*sigma3*log(bb1)
            ba13 = dcos(2.0*api*dd1)
              sf  = sqrt(ba)*ba13
            cellper = 960.0 + sf
c....................................................................
c     Create a representative cell in an asynchronous population
c....................................................................
            cd1 = duni()
              aphi = 1.0 - (log(2.0 - cd1)/log(2.0))
            asy_pop = aphi*cellper
               diff_time=cellper-asy_pop
             if(diff_time.le.aobs_time)then
               rest_time=aobs_time-diff_time
c                write(3002,*)rest_time
                else
c                  write(3001,*)cellper,asy_pop,'a'
                   endif
c....................................................................
                 write(*,*)cellper,asy_pop
c..................................................................
c                  Period define                 
c..................................................................
            G1_time=0.25*cellper
            S_time=0.75*cellper
            G2_time=0.85*cellper
        if(diff_time.gt.aobs_time) then    
          if(asy_pop.le.G1_time) then
             ak111=ak1*2
             ak41=ak4*2
             ak251=ak25*2
             ak291=ak29*2
             ak231=ak23*2
             ak271=ak27*2
             ak62=ak61*1
             ak102=ak101*1
c       write(2100,*) cellper,rest_time,'a', 'G1', ak62, ak102 
c             write(2200,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
            else
        if((asy_pop.gt.G1_time) .and. (asy_pop.le.S_time)) then
               ak111=ak1*2
               ak41=ak4*2
               ak231=ak23*2
               ak251=ak25*2
               ak271=ak27*2
               ak291=ak29*2
               ak62=ak61*2
               ak102=ak101*2
c       write(2100,*) cellper,rest_time,'a', 'S', ak62, ak102 
c             write(2300,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
              else
        if((asy_pop.gt.S_time) .and. (asy_pop.le.G2_time)) then
                ak111=ak1
                ak231=ak23
                ak41=ak4
                ak291=ak29
                ak251=ak25
                ak271=ak27
                ak62=ak61*2
                ak102=ak101*2
c       write(2100,*) cellper,rest_time,'a', 'G2', ak62, ak102 
c             write(2400,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
                   else
                    ak111=ak1
                    ak41=ak4
                    ak231=ak23
                    ak251=ak25
                    ak271=ak27
                    ak291=ak29
                    ak62=ak61*2
                    ak102=ak101*2
c      write(2100,*) cellper,rest_time,'a', 'M', ak62, ak102 
c             write(2500,*) cellper,rest_time, ak62, ak102 
c             write(3000,*) cellper,rest_time,ak62, ak102 
c               write(10,*) ak111,ak41,ak231,ak251,ak271,ak291     
               endif     
           endif
          endif 
        endif 
c......................................................................                  
c.....................................................................
c       How long to simulate a particular single cell?
c.....................................................................
             if(diff_time.gt.aobs_time)then
             ti=0.0
              t_f=aobs_time+1.0
               end_time=aobs_time
c.....................................................................
c      Initial conditions
c......................................................................
       Ga=0.0
       GaO=0.0
       GaS=0.0
       yM1=0.0
       yM2=16.0
       xP1=18.0
       zP2=33.0
       DP1=0.0
       DP2=1.0
       bOS=9.0
       ymO=35.0
       ymS=47.0
       bOct4=686.0
       Sox2=227.0
       ymiR=503.0
       ymiM1=9.0
       ymiM2=0.0
c                write(*,*)'ti=',ti,'t_f=',t_f
c.................................................................
              else
                ti=0.0
                t_f=diff_time+1.0
                  end_time=diff_time
       Ga=0.0
       GaO=0.0
       GaS=0.0
       yM1=0.0
       yM2=16.0
       xP1=18.0
       zP2=33.0
       DP1=0.0
       DP2=1.0
       bOS=9.0
       ymO=35.0
       ymS=47.0
       bOct4=686.0
       Sox2=227.0
       ymiR=503.0
       ymiM1=9.0
       ymiM2=0.0
c               write(*,*)'ti=',ti,'t_f=',t_f
              endif
c
              endif
c..................................................................
c
      nrxn=0
      faltucount=0
c
c     -----------------------------------------------------------------
c     STEP 1: CALCULATION OF a_nu=h_nu*c_nu AND a0
c     -----------------------------------------------------------------
 10     a_nu(1)=(((ak111/(ak2312**2+bOS**2))*(ak3-Ga)*ak90*bOS)*ak60)
 		a_nu(2)=((ak2-(ak393*ak382)+(ak394*ak382**2))*Ga*ak60)
 		a_nu(3)=(ak41*ak90*ak5*(ak3-Ga)*ak60)
 		a_nu(4)=(ak62*Ga*ak611*ak60)
 		a_nu(5)=(ak7*yM1*ak60)
 		a_nu(6)=(ak8*ymiR*yM1*ak60)
 		a_nu(7)=(ak9*ymiM1*ak60)
 		a_nu(8)=(ak15*ymiM1*ak60)
 		a_nu(9)=(ak102*ak1011*Ga*ak60)
 		a_nu(10)=(ak11*yM2*ak60)
 		a_nu(11)=(ak12*ymiR*yM2*ak60)
 		a_nu(12)=(ak13*ymiM2*ak60)
 		a_nu(13)=(ak15*ymiM2*ak60)
 		a_nu(14)=((ak14/(1+(ak390*ak381)+(ak395*ak382)))*ak60)
 		a_nu(15)=(ak15*yMiR*ak60)
 		a_nu(16)=(ak7*ymiM1*ak60)
 		a_nu(17)=(ak11*ymiM2*ak60)
 		a_nu(18)=(ak81*ymiM1*ak60)
 		a_nu(19)=(ak82*ymiM2*ak60)
 		a_nu(20)=(ak16*yM1*ak60)
 		a_nu(21)=(ak17*yM2*ak60)
 		a_nu(22)=(ak18*ymiM1*ak60)
 		a_nu(23)=(ak19*ymiM2*ak60)
 		a_nu(24)=(ak22*xP1*ak60)
 		a_nu(25)=(ak22*zP2*ak60)
 		a_nu(26)=(ak20*xP1*(xP1-1)*ak60)
 		a_nu(27)=(ak20*zP2*(zP2-1)*ak60)
 		a_nu(28)=(ak21*DP1*ak60)
 		a_nu(29)=(ak21*DP2*ak60)
 		a_nu(30)=(ak22*DP1*ak60)
 		a_nu(31)=(ak22*DP2*ak60)
 		a_nu(32)=(ak231*DP1*(ak26-GaO)*ak60)
 		a_nu(33)=(ak231*DP2*(ak26-GaO)*ak60)
 		a_nu(34)=(ak24*GaO*ak60)
 		a_nu(35)=(ak251*(ak26-GaO)*ak5*ak60)
 		a_nu(36)=(ak271*(ak30-GaS)*DP1*ak60)
 		a_nu(37)=(ak271*(ak30-GaS)*DP2*ak60)
 		a_nu(38)=(ak28*GaS*ak60)
 		a_nu(39)=(ak291*(ak30-GaS)*ak5*ak60)
 		a_nu(40)=(ak31*GaO*ak60)
 		a_nu(41)=(ak32*ymO*ak60)
 		a_nu(42)=(ak33*GaS*ak60)
 		a_nu(43)=(ak34*ymS*ak60)
 		a_nu(44)=(ak35*ymO*ak60)
 		a_nu(45)=(ak36*bOct4*ak60)
 		a_nu(46)=(ak39*bOct4*Sox2*ak60)
 		a_nu(47)=(ak40*bOS*ak60)
 		a_nu(48)=(ak37*ymS*ak60)
 		a_nu(49)=(ak38*Sox2*ak60)
 		a_nu(50)=(ak38*bOS*ak60)
 		a_nu(51)=(ak36*bOS*ak60) 		
 		a_nu(52)=(ak2311*Ga*DP1*ak60) 		
 		a_nu(53)=(ak2311*Ga*DP2*ak60)
c.........................................................................       
      sum1=0.0
       do 85 i=1,nu
         sum1=sum1+a_nu(i)
 85    continue
       a0=sum1
c       if(a0.eq.0.0)go to 75
c      write(*,*)'a0 =',a0

c     -----------------------------------------------------------------
c     STEP 2: CALCULATION OF tau AND mu
c     -----------------------------------------------------------------

c     tau is the infinitesimal time (continuous) during which a 
c     particular reaction takes place
c     mu determines the specific reaction channel

c      r1=rand(idum)
      r1=duni()
       if(r1.eq.0.0)then
         r1=0.001
         endif
c      call random_number(harvest=r1)
      tau=log(1.0/r1)/a0

c      r2=rand(idum)
       r2=duni()
c      call random_number(harvest=r2)
      r2a0=r2*a0

      sum2=0.0
      do 30 j=1,nu
         mu=j
         sum2=sum2+a_nu(j)
c	write(*,*)mu
         if(sum2.ge.r2a0)go to 100
 30   continue

c     -----------------------------------------------------------------
c     STEP 3: SAMPLING
c     -----------------------------------------------------------------
c......................................................................
c     this step updates the population of each reacting species
c......................................................................
 100  go to(101,102,103,104,105,106,107,108,109,110,111,112,113,
     $114,115,116,117,118,119,120,121,122,123,124,125,126,127,
     $128,129,130,131,132,133,134,135,136,137,138,139,140,141,
     $142,143,144,145,146,147,148,149,150,151,152,153),mu
 101     Ga=Ga+1.0
        go to 50
 102     Ga=Ga-1.0
        go to 50
 103     Ga=Ga+1.0
        go to 50
 104     yM1=yM1+1.0
        go to 50
 105     yM1=yM1-1.0
        go to 50
 106     ymiM1=ymiM1+1.0
 		 ymiR=ymiR-1.0
 		 yM1=yM1-1.0
        go to 50
 107     yM1=yM1+1.0
 		 ymiR=ymiR+1.0
        ymiM1=ymiM1-1.0
        go to 50
 108     ymiM1=ymiM1-1.0
 		 yM1=yM1+1.0
        go to 50
 109    yM2=yM2+1.0
        go to 50  
 110    yM2=yM2-1.0
        go to 50
 111    yM2=yM2-1.0
 		ymiR=ymiR-1.0
 		ymiM2=ymiM2+1.0
        go to 50
 112    yM2=yM2+1.0
 		ymiR=ymiR+1.0
 		ymiM2=ymiM2-1.0
        go to 50
 113     ymiM2=ymiM2-1.0
 		 yM2=yM2+1.0
        go to 50
 114    ymiR=ymiR+1.0
        go to 50     
 115    ymiR=ymiR-1.0
        go to 50
 116    ymiR=ymiR+1.0
 		ymiM1=ymiM1-1.0
        go to 50
 117    ymiR=ymiR+1.0
 		ymiM2=ymiM2-1.0
        go to 50       
 118     ymiM1=ymiM1-1.0
 		 ymiR=ymiR+1.0
        go to 50
 119    ymiM2=ymiM2-1.0
  		 ymiR=ymiR+1.0
        go to 50
 120    xP1=xP1+1.0
        go to 50     
 121    zP2=zP2+1.0
        go to 50
 122    xP1=xP1+1.0
        go to 50
 123    zP2=zP2+1.0
        go to 50
 124    xP1=xP1-1.0
        go to 50
 125    zP2=zP2-1.0
        go to 50
 126    xP1=xP1-2.0
        DP1=DP1+1.0
	    go to 50
 127    zP2=zP2-2.0
 		DP2=DP2+1.0
        go to 50
 128	DP1=DP1-1.0
 		xP1=xP1+2.0
 		go to 50
 129	DP2=DP2-1.0
 		zP2=zP2+2.0
 		go to 50
 130	DP1=DP1-1.0
 		xP1=xP1+1.0
 		go to 50
 131	DP2=DP2-1.0
 		zP2=zP2+1.0
 		go to 50  
 132    GaO=GaO+1.0
        go to 50
 133    GaO=GaO+1.0
        go to 50
 134    GaO=GaO-1.0
        go to 50
 135    GaO=GaO+1.0
        go to 50 
 136    GaS=GaS+1.0
        go to 50
 137    GaS=GaS+1.0
        go to 50
 138    GaS=GaS-1.0
	    go to 50
 139    GaS=GaS+1.0
        go to 50 
 140    ymO=ymO+1.0
        go to 50 
 141    ymO=ymO-1.0
        go to 50
 142    ymS=ymS+1.0
        go to 50
 143    ymS=ymS-1.0
        go to 50 
 144    bOct4=bOct4+1.0
        go to 50 
 145    bOct4=bOct4-1.0
        go to 50
 146    bOct4=bOct4-1.0
 		Sox2=Sox2-1.0
 		bOS=bOS+1.0
        go to 50
 147    bOS=bOS-1.0
 		bOct4=bOct4+1.0
 		Sox2=Sox2+1.0
        go to 50 
 148    Sox2=Sox2+1.0
        go to 50 
 149    Sox2=Sox2-1.0
        go to 50
 150    bOS=bOS-1.0
 		bOct4=bOct4+1.0
        go to 50
 151    Sox2=Sox2+1.0
 		bOS=bOS-1.0
        go to 50                
 152     Ga=Ga-1.0
        go to 50
 153     Ga=Ga-1.0
        go to 50
c---------------------------------------------------------------
 50      nrxn=nrxn+1
           ti=ti+tau
c          if(ti.ge.5.0)then
c           		agrof=3.0
c          endif
           
           
         if(nrxn.eq.1500000000)then
         nrxn=2
            faltucount=faltucount+1
c          write(*,*)faltucount,ti
         endif
c..............................................................
         if(ti.ge.end_time)then
c..............................................................
           if(ccount.eq.0.0)then
c..............................................................
          if(diff_time.ge.aobs_time)then
c..............................................................
              if(count.eq.0.0)then
                 count=count+1.0
                 X=(xP1+zP2+2*DP1+2*DP2)
                 yM=(yM1+yM2+ymiM1+ymiM2)
		    yMt1=yM1+ymiM1
		    yMt2=yM2+ymiM2
		    xPt1=xP1+2*DP1
		    zPt2=zP2+2*DP2 
		    write(2501,*)ti,xPt1
		    write(2502,*)ti,zPt2
		    write(2503,*)ti,yM1
		    write(2508,*)ti,yM2
		    write(2504,*)ti,bOS
		    write(2505,*)ti,bOct4
		    write(2506,*)ti,Sox2
		    write(2507,*)ti,X
		    write(2511,*)ti,yM
		    write(2513,*)ti,ymiR
		    write(2514,*)ti,ymiM1
		    write(2515,*)ti,ymiM2
		    x11(jj)=yM1+ymiM1
		    y11(jj)=yM2+ymiM2
		   x12(jj)=xP1+2*DP1
		   y12(jj)=zP2+2*DP2
		   x13(jj)=yM
		   y13(jj)=X	
c	 write(5003,*)x11(jj),y11(jj),x12(jj),y12(jj),x13(jj),y13(jj)	
		    write(2509,*)yM1,yM2,xP1,zP2,DP1,DP2
		    write(2512,*)yM,X,yMt1,yMt2,xPt1,zPt2
		    write(2516,*)yM1,ymiM1,yMt1,yM2,ymiM2,yMt2
		    write(2517,*)xP1,DP1,xPt1,zP2,DP2,zPt2		    		    		    		    
                 
c                    write(*,*)'ti=',ti,'t_f=',t_f
                   endif
c............................................................
                else
c.............................................................
                if(count.eq.0.0)then
               count=count+1.0
                 ixP11=xP1/2.0
                  xP11=ixP11
                  izP21=zP2/2.0
                   zP21=izP21
			 iDP11=DP1/2.0
			  DP11=iDP11
			   iDP21=DP2/2.0
			    DP21=iDP21
                   iyM51=yM1/2.0
                    yM51=iyM51
                    iyM52=yM2/2.0
                    yM52=iyM52
                     ibOS1=bOS/2.0
                       bOS1=ibOS1
                      iymO1=ymO/2.0
                        ymO1=iymO1
                       iymS1=ymS/2.0
                         ymS1=iymS1
                        ibOct41=bOct4/2.0
                         bOct41=ibOct41
                          iSox21=Sox2/2.0
                        Sox21=iSox21
                        	iymiR1=ymiR/2.0
                        	  ymiR1=iymiR1
                        	iymiM11=ymiM1/2.0
                        	 ymiM11=iymiM11
                       iymiM21=ymiM2/2.0
                         ymiM21=iymiM21	   
c       write(3003,*)xP11,zP21,DP11,DP21,yM51,yM52,bOS1
c                 write(2509,*)yM51,yM52
c        write(3004,*)ymO1,ymS1,bOct41,Sox21
c        write(3007,*)ymiR1,ymiM11,ymiM21
                     xP12=xP1-xP11
                      zP22=zP2-zP21
			     DP12=DP1-DP11
				DP22=DP2-DP21
                     yM61=yM1-yM51
                     yM62=yM2-yM52
                      bOS2=bOS-bOS1
                     ymO2=ymO-ymO1
                      ymS2=ymS-ymS1
                     bOct42=bOct4-bOct41
                      Sox22=Sox2-Sox21
                      ymiR2=ymiR-ymiR1
                      ymiM12=ymiM1-ymiM11
                      ymiM22=ymiM2-ymiM21
c              write(3005,*)xP12,zP22,DP12,DP22,yM61,yM62,bOS2
c              write(2509,*)yM61,yM62
c               write(3006,*)ymO2,ymS2,bOct42,Sox22
c        write(3008,*)ymiR2,ymiM12,ymiM22
               ccount=ccount+1.0

c                write(*,*)'ti=',ti,'t_f=',t_f
              endif
c...........................................................
                endif
c...........................................................
              else
c...........................................................
               if(bcount.eq.1.0)then
c...........................................................
                if(ccount.eq.1.0)then
		    if(count.eq.0.0)then
		    count=count+1.0
                 X=(xP1+zP2+2*DP1+2*DP2)
                 yM=(yM1+yM2+ymiM1+ymiM2)
		    yMt1=yM1+ymiM1
		    yMt2=yM2+ymiM2
		    xPt1=xP1+2*DP1
		    zPt2=zP2+2*DP2 
		    write(2501,*)ti,xPt1
		    write(2502,*)ti,zPt2
		    write(2503,*)ti,yM1
		    write(2508,*)ti,yM2
		    write(2504,*)ti,bOS
		    write(2505,*)ti,bOct4
		    write(2506,*)ti,Sox2
		    write(2507,*)ti,X
		    write(2511,*)ti,yM
		    write(2513,*)ti,ymiR
		    write(2514,*)ti,ymiM1
		    write(2515,*)ti,ymiM2
		    x11(jj)=yM1+ymiM1
		    y11(jj)=yM2+ymiM2
		   x12(jj)=xP1+2*DP1
		   y12(jj)=zP2+2*DP2
		   x13(jj)=yM
		   y13(jj)=X	
c	  write(5003,*)x11(jj),y11(jj),x12(jj),y12(jj),x13(jj),y13(jj)	
c		    write(*,*)x11(jj),y11(jj)
		    write(2510,*)yM1,yM2,xP1,zP2,DP1,DP2
		    write(2512,*)yM,X,yMt1,yMt2,xPt1,zPt2 
		    write(2516,*)yM1,ymiM1,yMt1,yM2,ymiM2,yMt2
		    write(2517,*)xP1,DP1,xPt1,zP2,DP2,zPt2		    		    		    
                    ccount=ccount+1.0
c                   write(*,*)'ti=',ti,'t_f=',t_f
                    endif
                      endif
c...........................................................
                  else
c...........................................................
                  if(bcount.eq.2.0)then
                    if(ccount.eq.2.0)then
                       if(count.eq.0.0)then
                     count=count+1.0
                 X=(xP1+zP2+2*DP1+2*DP2)
                 yM=(yM1+yM2+ymiM1+ymiM2)
		    yMt1=yM1+ymiM1
		    yMt2=yM2+ymiM2
		    xPt1=xP1+2*DP1
		    zPt2=zP2+2*DP2 
		    write(2501,*)ti,xPt1
		    write(2502,*)ti,zPt2
		    write(2503,*)ti,yM1
		    write(2508,*)ti,yM2
		    write(2504,*)ti,bOS
		    write(2505,*)ti,bOct4
		    write(2506,*)ti,Sox2
		    write(2507,*)ti,X
		    write(2511,*)ti,yM
		    write(2513,*)ti,ymiR
		    write(2514,*)ti,ymiM1
		    write(2515,*)ti,ymiM2
		    x11(jj)=yM1+ymiM1
		    y11(jj)=yM2+ymiM2
		   x12(jj)=xP1+2*DP1
		   y12(jj)=zP2+2*DP2
		   x13(jj)=yM
		   y13(jj)=X	
c	  write(5003,*)x11(jj),y11(jj),x12(jj),y12(jj),x13(jj),y13(jj)
		    write(2510,*)yM1,yM2,xP1,zP2,DP1,DP2
		    write(2512,*)yM,X,yMt1,yMt2,xPt1,zPt2
		    write(2516,*)yM1,ymiM1,yMt1,yM2,ymiM2,yMt2
		    write(2517,*)xP1,DP1,xPt1,zP2,DP2,zPt2		    		    		    
                               ccount=0.0
                              bcount=0.0
c                          write(*,*)'ti=',ti,'t_f=',t_f
                          endif
                        endif
                        endif
c...........................................................
                          endif
c...........................................................
                            endif
c.................................................................
                             endif
c..................................................................
c....................................................................
c       To draw stochastic profile at any instance
c......................................................................
        if(jj.eq.1)then
	its=mod(nrxn,1000)
	if(its.eq.0)then
c        if(ti .ge. 200000)then
c          tscaled=ti/10000.0
c         if(tscaled .ge. 20)then
c          ts=(tscaled-20)
                 X=(xP1+zP2+2*DP1+2*DP2)
                 yM=(yM1+yM2+ymiM1+ymiM2)
		    yMt1=yM1+ymiM1
		    yMt2=yM2+ymiM2
		    xPt1=xP1+2*DP1
		    zPt2=zP2+2*DP2 
c        write(2401,*)ti,Ga
c        write(2402,*)ti,GaO
c        write(2403,*)ti,GaS
c        write(2404,*)ti,xPt1
c        write(2405,*)ti,zPt2
c        write(2406,*)ti,yM
c        write(2407,*)ti,bOS
c        write(2408,*)ti,bOct4
c        write(2409,*)ti,Sox2
c        write(2410,*)ti,ymO
c        write(2411,*)ti,ymS
c        write(2412,*)ti,X
c        write(2413,*)ti,yMt1
c        write(2414,*)ti,yMt2
c        write(2415,*)ti,ymiR
c        write(2416,*)ti,ymiM1
c        write(2417,*)ti,ymiM2

c        endif
c         endif
	   endif
           endif
c.......................................................................
c         Time loop continues
c......................................................................
            if(ti.lt.t_f)then
			go to 10		
    		endif
    		If(ccount.eq.1.0)then
    		count=0.0
    		go to 355
    		endif
c           write(*,*)'ti=',ti,'t_f=',t_f     
 300    continue
c.......................................................................
c         Noise calculation 
c......................................................................
  
            sum13=0
            sum14=0
            sum3=0
            sum4=0
            sum5=0
            sum6=0
            sum7=0
            sum8=0
            sum9=0
            sum10=0
            sum11=0
            sum12=0    
        do 5 jl=1,nk  
c        read(2509,*) x11(ll),y11(ll)  
c      write(1000,*) x11(jl),y11(jl),x12(jl),y12(jl), x13(jl), y13(jl)
c.....................................................................
c        mRNA noise calculation 
c......................................................................
            sqx1(jl)= x11(jl)*x11(jl)
            sqy1(jl)= y11(jl)*y11(jl)
            multi1(jl)= x11(jl)*y11(jl)
            diff1(jl)= x11(jl)-y11(jl)
            sq1(jl)= diff1(jl)*diff1(jl)
            sqt1(jl)= sqx1(jl)+sqy1(jl)
            sum13=sum13+x11(jl)
            sum14=sum14+y11(jl)
            sum3=sum3+multi1(jl)
            sum4=sum4+sq1(jl)
            sum5=sum5+sqt1(jl)
            sum11=sum11+x13(jl)
c.....................................................................
c        protein noise calculation 
c.....................................................................
		sqx2(jl)=x12(jl)*x12(jl)
		sqy2(jl)=y12(jl)*y12(jl)
            multi2(jl)= x12(jl)*y12(jl)
            diff2(jl)= x12(jl)-y12(jl)
            sq2(jl)= diff2(jl)*diff2(jl)
            sqt2(jl)= sqx2(jl)+sqy2(jl)
            sum6=sum6+x12(jl)
            sum7=sum7+y12(jl)
            sum8=sum8+multi2(jl)
            sum9=sum9+sq2(jl)
            sum10=sum10+sqt2(jl)
            sum12=sum12+y13(jl)
c......................................................................
    5       continue
c......................................................................
c		mRNA noise
c......................................................................
            A1=sum13/nk
            B1=sum14/nk
            C1=sum4/nk
            D1=sum3/nk
            F1=sum5/nk
            G1=2*A1*B1
            H1=A1*B1
            Q1=(C1/G1)
            S1=(D1-H1)/H1
            T1=(F1-G1)/G1
            W1=Sum11/nk
c            write(*,*) A1, B1, C1, D1, F1, G1, H1
            write(*,*)Q1,S1,T1, W1
            write(5001,*)Q1,S1,T1,W1
c......................................................................
c		protein noise
c......................................................................
            A2=sum6/nk
            B2=sum7/nk
            C2=sum9/nk
            D2=sum8/nk
            F2=sum10/nk
            G2=2*A2*B2
            H2=A2*B2
            Q2=(C2/G2)
            S2=(D2-H2)/H2
            T2=(F2-G2)/G2
            W2=sum12/nk
            write(*,*)Q2,S2,T2, W2
            write(5002,*)Q2,S2,T2, W2
c........................................................................
       stop
       end
c........................................................................
c      Subroutine for random no generator by Lane Watson
c.........................................................................
            DOUBLE PRECISION FUNCTION DUNI()
c            FUNCTION DUNI()
C***BEGIN PROLOGUE  DUNI
C***DATE WRITTEN   880714 (YYMMDD)
C***REVISION DATE  880714 (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  THIS ROUTINE GENERATES DOUBLE PRECISION UNIFORM
C             RANDOM NUMBERS ON [0,1)
C***DESCRIPTION
C        COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C
C       USAGE:
C              TO INITIALIZE THE GENERATOR
C                   USEED = DUSTAR(ISEED)
C               WHERE: ISEED IS ANY NONZERO INTEGER
C                  WILL RETURN FLOATING POINT VALUE OF ISEED.
C
C               SUBSEQUENTLY
C                       U = DUNI()
C                  WILL RETURN A REAL UNIFORM ON [0,1)
C
C                ONE INITIALIZATION IS NECESSARY, BUT ANY NUMBER OF EVALUATIONS
C                  OF DUNI IN ANY ORDER, ARE ALLOWED.
C
C           NOTE: DEPENDING UPON THE VALUE OF K (SEE BELOW), THE OUTPUT
C                       OF DUNI MAY DIFFER FROM ONE MACHINE TO ANOTHER.
C
C           TYPICAL USAGE:
C
C               DOUBLE PRECISION U,DUNI,DUSTAR,USEED
C               INTEGER ISEED
CC                 SET SEED
C               ISEED = 305
C               USEED = DUSTAR(ISEED)
C               DO 1 I = 1,1000
C                   U = DUNI()
C             1 CONTINUE
CC                 NOTE: IF K=47 (THE DEFAULT, SEE BELOW) THE OUTPUT VALUE OF
CC                           U WILL BE 0.812053811384E-01...
C               WRITE(*,*) U
C               END
C
C          NOTE ON PORTABILITY: USERS CAN CHOOSE TO RUN DUNI IN ITS DEFAULT
C               MODE (REQUIRING NO USER ACTION) WHICH WILL GENERATE THE SAME
C               SEQUENCE OF NUMBERS ON ANY COMPUTER SUPPORTING FLOATING POINT
C               NUMBERS WITH AT LEAST 47 BIT MANTISSAS, OR IN A MODE THAT
C               WILL GENERATE NUMBERS WITH A LONGER PERIOD ON COMPUTERS WITH
C               LARGER MANTISSAS.
C          TO EXERCISE THIS OPTION:  B E F O R E  INVOKING DUSTAR INSERT
C               THE INSTRUCTION        UBITS = DUNIB(K)      K >= 47
C               WHERE K IS THE NUMBER OF BITS IN THE MANTISSA OF YOUR FLOATING

C               POINT WORD (K=96 FOR CRAY, CYBER 205). DUNIB RETURNS THE
C               FLOATING POINT VALUE OF K THAT IT ACTUALLY USED.
C                    K INPUT AS .LE. 47, THEN UBITS=47.
C                    K INPUT AS .GT. 47, THEN UBITS=FLOAT(K)
C               IF K>47 THE SEQUENCE OF NUMBERS GENERATED BY DUNI MAY DIFFER
C               FROM ONE COMPUTER TO ANOTHER.
C
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE DUNI
      DOUBLE PRECISION CSAVE,CD,CM
      PARAMETER(
     *    CSAVE=0.9162596898123D+13/0.140737488355328D+15,
     *    CD=0.76543212345678D+14/0.140737488355328D+15,
     *    CM=0.140737488355213D+15/0.140737488355328D+15)
C                            2**47=0.140737488355328D+15
      DOUBLE PRECISION U(17),S,T,DUSTAR,C,DUNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED
C
      SAVE U,I,J,K,C
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH ISEED=305 AND K=96.
      DATA U/
     *0.471960981577884755837789724978D+00,
     *0.930323453205669578433639632431D+00,
     *0.110161790933730836587127944899D+00,
     *0.571501996273139518362638757010D-01,
     *0.402467554779738266237538503137D+00,
     *0.451181953427459489458279456915D+00,
     *0.296076152342721102174129954053D+00,
     *0.128202189325888116466879622359D-01,
     *0.314274693850973603980853259266D+00,
     *0.335521366752294932468163594171D-02,
     *0.488685045200439371607850367840D+00,
     *0.195470426865656758693860613516D+00,
     *0.864162706791773556901599326053D+00,
     *0.335505955815259203596381170316D+00,
     *0.377190200199058085469526470541D+00,
     *0.400780392114818314671676525916D+00,
     *0.374224214182207466262750307281D+00/
      DATA I,J,K,C/17,5,47,CSAVE/
C
C   BASIC GENERATOR IS FIBONACCI
C
      DUNI = U(I)-U(J)
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      U(I) = DUNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
C   SECOND GENERATOR IS CONGRUENTIAL
C
      C = C-CD
      IF(C.LT.0.0D0) C=C+CM
C
C   COMBINATION GENERATOR
C
      DUNI = DUNI-C
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      RETURN
C
      ENTRY DUSTAR(ISEED)
C
C          SET UP ...
C          CONVERT ISEED TO FOUR SMALLISH POSITIVE INTEGERS.
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED.
C
        DO 2 II = 1,17
          S = 0.0D0
          T = 0.5D0
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER K BITS, WHERE K IS DEFAULTED TO 47 BUT CAN
C               BE CHANGED BY USER CALL TO DUNIB(K)
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179)
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169)
                  IF(MOD(L1*M1,64).GE.32)S=S+T
    3             T = 0.5D0*T
    2   U(II) = S
        DUSTAR = FLOAT(ISEED)
        RETURN
C
      ENTRY DUNIB(KK)
        IF(KK.LE.47)THEN
             K=47
        ELSE
             K=KK
        ENDIF
        DUNIB=FLOAT(K)
      END

