module constantdata
	real*8,parameter::pi=dacos(-1.d0)
end


program main
	use constantdata
	implicit none
	integer i,j,l,r,ml,sl,mr,sr
	!i,j,l,r为高斯波函数循环变量;
	integer x(4),p
    real*8 symbol,ss(6,5)
	!x(4)用于读取右矢波函数粒子顺序循环变量；p是右矢的个数,考虑到12对称,34对称所以只有6个；symbol是符号；
	!ss(6,5)是用于存放反对称化后的6个波函数和symbol
	
	real*8 s(1000),si0,b0,miu(1000),b(1000)
	!12-34相对运动的高斯波函数中心个数*12团内的宽度个数=高斯个数
	real*8 N(1000,1000),H(1000,1000),H1(1000,1000),H2(1000,1000),Nbar(1000,1000),NbarT(1000,1000),Hbar(1000,1000)
	!N矩阵元,H矩阵元=H1(势能)+H2(动能)+4*mass*N
	real*8 vv(1000)
	!存放插值后的势能函数vv(r)
	integer info
	real*8 mass,hba
	!miu是高斯波函数宽度,mass是中子质量,hba是普朗克常数
	real*8 w(1000),work(1000)
    real*8 functionn,v12,v13,v14,v23,v24,v34,kinetic
	integer gaosinum,gaosinummax,sgaosinum,miugaosinum,gaosinumbar
	character(len=3) fname
	!filename用来产生文件名
	integer nbari
    real*8  Wmin
    !用于读取人为设定的最小本征能量

	call chazhi(vv)

	
	open(10,file='s.dat',status='old')
        !read(10,*) miu

		read(10,*) Wmin
        read(10,*) mass
		read(10,*) hba
		read(10,*) ((ss(i,j),j=1,5),i=1,6)
		print*, "Please input gaosinummax;(must be greater than 1)"
		read(*,*) gaosinummax
		!gaosinum 要求运行程序的时候输入
	close(10)

do sgaosinum=1,gaosinummax
do miugaosinum=2,gaosinummax
	gaosinum=sgaosinum*miugaosinum
		b0=90.d0/(miugaosinum+1)
        si0=20.d0/(sgaosinum+1)
		!中子的渐近区在无穷远;
		!有效波函数为0~3;但是为了排除人为选择对结果的影响,将波函数中心位置范围扩大到6;
    	N=0.d0
	    H=0.d0
	    H1=0.d0
	    H2=0.d0
		do i=1,miugaosinum
			b(i)=b0*i
			miu(i)=0.5d0/(b(i)**2.d0)
		enddo
		do i=1,sgaosinum
			s(i)=si0*i
		enddo
		!获得sgaosinum个团间高斯中心;miugaosinum个团内高斯宽度;
		!可以通过加个常数实现所有高斯宽度的,高斯中心的整体平移;
		
	do p=1,6
		symbol=ss(p,1)
		do i=1,4
		x(i)=ss(p,1+i)
		!读入右矢波函数,及其符号
		enddo
		
		i=0
		do ml=1,miugaosinum
		do sl=1,sgaosinum
			i=i+1
			j=0
				do mr=1,miugaosinum
				do sr=1,sgaosinum
					j=j+1
					N(i,j)=N(i,j)+symbol*(4.d0*miu(ml)/Pi)**3.d0*(4.d0*miu(mr)/Pi)**3.d0*functionN(x,s(sl),s(sr),miu(ml),miu(mr))
					H1(i,j)=H1(i,j)+symbol*(4.d0*miu(ml)/Pi)**3.d0*(4.d0*miu(mr)/Pi)**3.d0*(V12(x,s(sl),s(sr),miu(ml),miu(mr),vv)+V13(x,s(sl),s(sr),miu(ml),miu(mr),vv)+V14(x,s(sl),s(sr),miu(ml),miu(mr),vv)+V23(x,s(sl),s(sr),miu(ml),miu(mr),vv)+V24(x,s(sl),s(sr),miu(ml),miu(mr),vv)+V34(x,s(sl),s(sr),miu(ml),miu(mr),vv))
					H2(i,j)=H2(i,j)-hba**2.d0*symbol*(4.d0*miu(ml)/Pi)**3.d0*(4.d0*miu(mr)/Pi)**3.d0/(8.d0*mass)*kinetic(x,s(sl),s(sr),miu(ml),miu(mr))
				enddo	
				enddo
		enddo
		enddo
	enddo
	H=H1+H2+4.d0*mass*N
	
	write(fname,'(i3)') gaosinum
	open(4,file='N'//fname//'.dat') 
	!open(41,file='H1'//fname//'.dat')
	!open(42,file='H2'//fname//'.dat')
	!open(43,file='H'//fname//'.dat')
		!output N matrix
		!output H1 matrix
		!output H2 matrix
		!output H matrix
		do i=1,gaosinum
			do j=1,gaosinum
				write(4,"(f30.15,\)") N(i,j)
	!			write(41,"(f30.15,\)") H1(i,j)
	!			write(42,"(f30.15,\)") H2(i,j)
	!			write(43,"(f30.15,\)") H(i,j)
			enddo
			write(4,*)
	!		write(41,*)
	!		write(42,*)
	!		write(43,*)
		enddo
	close(4)
	!close(41)
	!close(42)
	!close(43)
	
    w=0.d0
    call dsyev('V','U',gaosinum,N,1000,w,work,3000,info)
	open(5,file='N_eigenvalue_'//fname//'.dat')
        !output N eigenvalue
		do i=1,gaosinum
			write(5,*) w(i)
		enddo
	close(5)
	
	if(w(1).lt.(10.d0**(-12.d0)))then
	!去掉本征值小于10**(-12)的本征矢,然后再重新构造H;再解本征能量;
		Nbar=0.d0
		Hbar=0.d0
		
		gaosinumbar=0
		do i=1,gaosinum
			if(w(i).gt.(10.d0**(-12.d0)))then
				gaosinumbar=gaosinumbar+1
				do j=1,gaosinum
					Nbar(gaosinumbar,j)=N(i,j)
					NbarT(j,gaosinumbar)=N(i,j)
				enddo
			endif
		enddo
		Hbar=matmul(matmul(NbarT,H),Nbar)
		w=0.d0
		call dsygv(1,'V','U',gaosinumbar,Hbar,1000,Nbar,1000,w,work,3000,info)
		
		open(6,file='en.dat')
        !output energy level
		!if(w(1).lt.Wmin)then
		!如果解出比参考能量小的能量本征值,输出;
			write(6,*) gaosinum
			do i=1,miugaosinum
			    write(6,"(f30.15,\)") miu(i)
            enddo
            do i=1,sgaosinum
			    write(6,"(f30.15,\)") s(i)
            enddo
			do i=1,gaosinum
			write(6,*)
				if(w(i).gt.0.d0.and.w(i).lt.10.d0**4.d0)then
					write(6,"(f30.15,\)") w(i)
				endif
			enddo
				write(6,*)
		!endif
	
	else
		open(50,file='N'//fname//'.dat',status='old')
			read(50,*) ((N(i,j),j=1,gaosinum),i=1,gaosinum)
		close(50)
	
		w=0.d0
		call dsygv(1,'V','U',gaosinum,H,1000,N,1000,w,work,3000,info)
		
		open(6,file='en.dat')
        !output energy level
		!if(w(1).lt.Wmin)then
		!如果解出比参考能量小的能量本征值,输出;
			write(6,*) gaosinum
            do i=1,miugaosinum
			    write(6,"(f30.15,\)") miu(i)
            enddo
            do i=1,sgaosinum
			    write(6,"(f30.15,\)") s(i)
            enddo
			write(6,*)
			do i=1,gaosinum
				if(w(i).gt.0.d0.and.w(i).lt.10.d0**4.d0)then
					write(6,"(f30.15,\)") w(i)
				endif
			enddo
				write(6,*)
		!endif
	endif
enddo
enddo	
	close(6)	
end

real*8 function functionN(x,sl,sr,miul,miur)
	use constantdata
    implicit none
    integer x(4),i
	real*8 miul,miur,sl,sr
	functionN=1.d0
	!四个粒子的overlap可对应“单独”计算
	do i=1,4
		!11 22 33 44 12 21 34 43
		if ((i.eq.x(i)).or.(i.eq.1.and.x(i).eq.2).or.(i.eq.2.and.x(i).eq.1).or.(i.eq.3.and.x(i).eq.4).or.(i.eq.4.and.x(i).eq.3)) then
			functionN=functionN*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		else
			!13 14 31 41 23 24 32 42
			functionN=functionN*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	enddo

end

real*8 function V12(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	v12=0.d0
	r=0.d0
	dr=0.003d0

		!12 21
		if ((x(1).eq.1.and.x(2).eq.2).or.(x(1).eq.2.and.x(2).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V12=V12+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*dr/6.d0
				else
					V12=V12+4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))
			enddo
			v12=v12*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl-sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl-sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(1).eq.1.and.x(2).eq.3).or.(x(1).eq.1.and.x(2).eq.4).or.(x(1).eq.2.and.x(2).eq.3).or.(x(1).eq.2.and.x(2).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V12=V12+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)*dr/6.d0
				else
					V12=V12+(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)
			enddo
			v12=v12*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(1).eq.3.and.x(2).eq.1).or.(x(1).eq.4.and.x(2).eq.1).or.(x(1).eq.3.and.x(2).eq.2).or.(x(1).eq.4.and.x(2).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V12=V12+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)*dr/6.d0
				else
					V12=V12+(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)
			enddo
			v12=v12*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(1).eq.3.and.x(2).eq.4).or.(x(1).eq.4.and.x(2).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V12=V12+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*dr/6.d0
				else
					V12=V12+4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))
			enddo
			v12=v12*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl+sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl+sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	!overlap
	do i=3,4
	   !12
	   if ((i.eq.x(i)).or.(i.eq.3.and.x(i).eq.4).or.(i.eq.4.and.x(i).eq.3)) then
			V12=v12*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			v12=v12*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	enddo
end

real*8 function V13(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	V13=0.d0
	r=0.d0
	dr=0.003d0
	
		!12 21
		if ((x(1).eq.1.and.x(3).eq.2).or.(x(1).eq.2.and.x(3).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V13=V13+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V13=V13+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V13=V13*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(1).eq.1.and.x(3).eq.3).or.(x(1).eq.1.and.x(3).eq.4).or.(x(1).eq.2.and.x(3).eq.3).or.(x(1).eq.2.and.x(3).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V13=V13+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*dr/6.d0
				else
					V13=V13+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)
                enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V13=V13*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(1).eq.3.and.x(3).eq.1).or.(x(1).eq.4.and.x(3).eq.1).or.(x(1).eq.3.and.x(3).eq.2).or.(x(1).eq.4.and.x(3).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V13=V13+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*dr/6.d0
				else
					V13=V13+4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V13=V13*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(1).eq.3.and.x(3).eq.4).or.(x(1).eq.4.and.x(3).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V13=V13+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V13=V13+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V13=V13*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	i=2
	   !12
	   if ((x(i).eq.2).or.(x(i).eq.1)) then
			V13=V13*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V13=V13*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	i=4
	   !34
	   if ((x(i).eq.4).or.(x(i).eq.3)) then
			V13=V13*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V13=V13*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
end

real*8 function V14(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	V14=0.d0
	r=0.d0
	dr=0.003d0
	
		!12 21
		if ((x(1).eq.1.and.x(3).eq.2).or.(x(1).eq.2.and.x(3).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V14=V14+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V14=V14+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V14=V14*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(1).eq.1.and.x(3).eq.3).or.(x(1).eq.1.and.x(3).eq.4).or.(x(1).eq.2.and.x(3).eq.3).or.(x(1).eq.2.and.x(3).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V14=V14+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*dr/6.d0
				else
					V14=V14+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V14=V14*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(1).eq.3.and.x(3).eq.1).or.(x(1).eq.4.and.x(3).eq.1).or.(x(1).eq.3.and.x(3).eq.2).or.(x(1).eq.4.and.x(3).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V14=V14+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*dr/6.d0
				else
					V14=V14+4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V14=V14*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(1).eq.3.and.x(3).eq.4).or.(x(1).eq.4.and.x(3).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V14=V14+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V14=V14+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V14=V14*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	i=2
	   !12
	   if ((x(i).eq.2).or.(x(i).eq.1)) then
			V14=V14*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V14=V14*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	i=3
	   !34
	   if ((x(i).eq.4).or.(x(i).eq.3)) then
			V14=V14*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V14=V14*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
end

real*8 function V23(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	V23=0.d0
	r=0.d0
	dr=0.003d0
	
		!12 21
		if ((x(1).eq.1.and.x(3).eq.2).or.(x(1).eq.2.and.x(3).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V23=V23+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V23=V23+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V23=V23*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(1).eq.1.and.x(3).eq.3).or.(x(1).eq.1.and.x(3).eq.4).or.(x(1).eq.2.and.x(3).eq.3).or.(x(1).eq.2.and.x(3).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V23=V23+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*dr/6.d0
				else
					V23=V23+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V23=V23*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(1).eq.3.and.x(3).eq.1).or.(x(1).eq.4.and.x(3).eq.1).or.(x(1).eq.3.and.x(3).eq.2).or.(x(1).eq.4.and.x(3).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V23=V23+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*dr/6.d0
				else
					V23=V23+4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V23=V23*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(1).eq.3.and.x(3).eq.4).or.(x(1).eq.4.and.x(3).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V23=V23+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V23=V23+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V23=V23*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	i=1
	   !12
	   if ((x(i).eq.2).or.(x(i).eq.1)) then
			V23=V23*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V23=V23*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	i=4
	   !34
	   if ((x(i).eq.4).or.(x(i).eq.3)) then
			V23=V23*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V23=V23*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
end


real*8 function V24(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	V24=0.d0
	r=0.d0
	dr=0.003d0
	
		!12 21
		if ((x(1).eq.1.and.x(3).eq.2).or.(x(1).eq.2.and.x(3).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V24=V24+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V24=V24+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V24=V24*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(1).eq.1.and.x(3).eq.3).or.(x(1).eq.1.and.x(3).eq.4).or.(x(1).eq.2.and.x(3).eq.3).or.(x(1).eq.2.and.x(3).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V24=V24+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*dr/6.d0
				else
					V24=V24+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul-sr*miur)))*Exp(-miul*(-sl+r)**2.d0-miur*(-sr+r)**2.d0))/(-sl*miul-sr*miur)
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V24=V24*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(1).eq.3.and.x(3).eq.1).or.(x(1).eq.4.and.x(3).eq.1).or.(x(1).eq.3.and.x(3).eq.2).or.(x(1).eq.4.and.x(3).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V24=V24+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*dr/6.d0
				else
					V24=V24+4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-miul*(sl**2.d0+r**2.d0)-miur*(sr**2.d0+r**2.d0))
			enddo
			!eee=Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V24=V24*Exp(-4.d0*(miul*miur)/(miul+miur)*(0.d0-0.d0)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(1).eq.3.and.x(3).eq.4).or.(x(1).eq.4.and.x(3).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V24=V24+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*dr/6.d0
				else
					V24=V24+(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sl*miul)))*Exp(-miul*(-sl+r)**2.d0-miur*(r)**2.d0))/(-sl*miul)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V24=V24*Exp(-1.d0*(miul*miur)/(miul+miur)*(sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	i=1
	   !12
	   if ((x(i).eq.2).or.(x(i).eq.1)) then
			V24=V24*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V24=V24*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	i=3
	   !34
	   if ((x(i).eq.4).or.(x(i).eq.3)) then
			V24=V24*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			V24=V24*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
end


real*8 function v34(x,sl,sr,miul,miur,vv)
	use constantdata
    implicit none
    integer x(4),i
	real*8 r,dr,miul,miur,sl,sr,vv(1000)
	v34=0.d0
	r=0.d0
	dr=0.003d0
	
		!12 21
		if ((x(3).eq.1.and.x(4).eq.2).or.(x(3).eq.2.and.x(4).eq.1)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V34=V34+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*dr/6.d0
				else
					V34=V34+4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl+sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			v34=v34*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl+sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!13 14 23 24
		if ((x(3).eq.1.and.x(4).eq.3).or.(x(3).eq.1.and.x(4).eq.4).or.(x(3).eq.2.and.x(4).eq.3).or.(x(3).eq.2.and.x(4).eq.4)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V34=V34+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)*dr/6.d0
				else
					V34=V34+(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(-sr*miur)))*Exp(-miul*(r)**2.d0-miur*(-sr+r)**2.d0))/(-sr*miur)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V34=V34*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!31 41 32 42
		if ((x(3).eq.3.and.x(4).eq.1).or.(x(3).eq.4.and.x(4).eq.1).or.(x(3).eq.3.and.x(4).eq.2).or.(x(3).eq.4.and.x(4).eq.2)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V34=V34+(vv(i-1)+4.d0*vv(i)+vv(i+1))*(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)*dr/6.d0
				else
					V34=V34+(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)*vv(i)*dr
				endif
				!eee=(pi*r*(-1.d0+Exp(4.d0*r*(sr*miur)))*Exp(-miul*(r)**2.d0-miur*(sr+r)**2.d0))/(sr*miur)
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V34=V34*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif
		!34 43
		if ((x(3).eq.3.and.x(4).eq.4).or.(x(3).eq.4.and.x(4).eq.3)) then
			!辛普森积分
			do i=1,1000
				r=i*dr
				if(i.ne.1000.and.i.ne.1.and.r.gt.dr)then
					V34=V34+(vv(i-1)+4.d0*vv(i)+vv(i+1))*4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*dr/6.d0
				else
					V34=V34+4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))*vv(i)*dr
				endif
				!eee=4.d0*pi*r**2.d0*Exp(-1.d0*r**2.d0*(miul+miur))
			enddo
			!eee=Exp(-1.d0*(miul*miur)/(miul+miur)*(sl-sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
			V34=V34*Exp(-1.d0*(miul*miur)/(miul+miur)*(sl-sr)**2.d0)*(pi**1.5d0)/(8.d0*(miul+miur)**1.5d0)
		endif

	do i=1,2
	   if ((i.eq.x(i)).or.(i.eq.1.and.x(i).eq.2).or.(i.eq.2.and.x(i).eq.1)) then
			v34=v34*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
	   else
			v34=v34*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
		endif
	enddo
end


    
real*8 function kinetic(x,sl,sr,miul,miur)
	use constantdata
    implicit none
    integer x(4),i,j,t
	real*8 miul,miur,sl,sr,kinetic1,kinetic2,kinetic3,kinetic4,overlap
	!kinetic1是pi**2.d0;kinetic2*kinetic3是pi*pj;overlap是i,j之外的两个波函数的耦合；
	kinetic=0.d0
	kinetic1=0.d0
		do i=1,4
			overlap=1.d0
			if(i.eq.1)then
				do j=2,4
					!11 22 33 44 12 21 34 43
					if ((j.eq.x(j)).or.(j.eq.1.and.x(j).eq.2).or.(j.eq.2.and.x(j).eq.1).or.(j.eq.3.and.x(j).eq.4).or.(j.eq.4.and.x(j).eq.3)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
					!13 14 31 41 23 24 32 42
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				enddo
			else
				if(i.eq.4)then
					do j=1,3
						!11 22 33 44 12 21 34 43
						if ((j.eq.x(j)).or.(j.eq.1.and.x(j).eq.2).or.(j.eq.2.and.x(j).eq.1).or.(j.eq.3.and.x(j).eq.4).or.(j.eq.4.and.x(j).eq.3)) then
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						else
						!13 14 31 41 23 24 32 42
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						endif
					enddo
				else
					do j=1,i-1
						!11 22 33 44 12 21 34 43
						if ((j.eq.x(j)).or.(j.eq.1.and.x(j).eq.2).or.(j.eq.2.and.x(j).eq.1).or.(j.eq.3.and.x(j).eq.4).or.(j.eq.4.and.x(j).eq.3)) then
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						else
						!13 14 31 41 23 24 32 42
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						endif
					enddo
					do j=i+1,4
						!11 22 33 44 12 21 34 43
						if ((j.eq.x(j)).or.(j.eq.1.and.x(j).eq.2).or.(j.eq.2.and.x(j).eq.1).or.(j.eq.3.and.x(j).eq.4).or.(j.eq.4.and.x(j).eq.3)) then
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						else
						!13 14 31 41 23 24 32 42
							overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
						endif
					enddo
				endif
			endif
			
			if ((i.eq.x(i)).or.(i.eq.1.and.x(i).eq.2).or.(i.eq.2.and.x(i).eq.1).or.(i.eq.3.and.x(i).eq.4).or.(i.eq.4.and.x(i).eq.3)) then
				kinetic1=kinetic1+2.d0**0.5d0*pi**1.5d0/(miul+miur)**3.5d0*Exp(-2.d0*miul*miur/(miul+miur)*(0.d0-0.d0)**2.d0)*(-3.d0*miur+miul*(-3.d0+4.d0*miur*(0.d0-0.d0)**2.d0))*overlap
			else
				kinetic1=kinetic1+2.d0**0.5d0*pi**1.5d0/(miul+miur)**3.5d0*Exp(-2.d0*miul*miur/(miul+miur)*(sl-sr)**2.d0)*(-3.d0*miur+miul*(-3.d0+4.d0*miur*(sl-sr)**2.d0))*overlap
			endif
		enddo
				kinetic1=3.d0*kinetic1
				
				
		do i=1,3
		do j=i+1,4
            
			overlap=1.d0
			!12
			if(i.eq.1.and.j.eq.2) then
				do t=3,4
					if ((t.eq.x(t)).or.(t.eq.3.and.x(t).eq.4).or.(t.eq.4.and.x(t).eq.3)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				enddo
			endif
			!13
			if(i.eq.1.and.j.eq.3) then
				t=2
					if ((t.eq.x(t)).or.(x(t).eq.1)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				t=4
					if ((t.eq.x(t)).or.(x(t).eq.3)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
			endif
			!23
			if(i.eq.2.and.j.eq.3) then
				t=1
					if ((t.eq.x(t)).or.(x(t).eq.2)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				t=4
					if ((t.eq.x(t)).or.(x(t).eq.3)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
			endif
			!14 24 34
			if(i.eq.1.and.j.eq.4) then
				do t=2,3
					if ((t.eq.x(t)).or.(t.eq.2.and.x(t).eq.1).or.(t.eq.3.and.x(t).eq.4)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				enddo
			endif
			if(i.eq.2.and.j.eq.4) then
				t=1
					if ((t.eq.x(t)).or.(x(t).eq.2)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
				t=3
					if ((t.eq.x(t)).or.(x(t).eq.4)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
			endif
			if(i.eq.3.and.j.eq.4) then
					do t=1,2
						if ((t.eq.x(t)).or.(t.eq.2.and.x(t).eq.1).or.(t.eq.1.and.x(t).eq.2)) then
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl-sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					else
						overlap=overlap*exp(-0.25d0*(miul*miur/(miul+miur))*(sl+sr)**2.d0)*(Pi**1.5d0)/((2.d0*(miul+miur))**1.5d0)
					endif
					enddo
			endif
			
			
			if((i.eq.x(i)).or.(i.eq.1.and.x(i).eq.2).or.(i.eq.2.and.x(i).eq.1).or.(i.eq.3.and.x(i).eq.4).or.(i.eq.4.and.x(i).eq.3)) then
			!11 22 33 44 12 21 34 43
				kinetic2=2.d0**0.5d0*pi**1.5d0*miul*miur/(miul+miur)**2.5d0*(0.d0-0.d0)*Exp(-2.d0*miul*miur/(miul+miur)*(0.d0-0.d0))
			else
			!13 23 14 24 31 32 41 42
				kinetic2=2.d0**0.5d0*pi**1.5d0*miul*miur/(miul+miur)**2.5d0*(sl-sr)*Exp(-2.d0*miul*miur/(miul+miur)*(sl-sr))
			endif
			
			if((j.eq.x(j)).or.(j.eq.1.and.x(j).eq.2).or.(j.eq.2.and.x(j).eq.1).or.(j.eq.3.and.x(j).eq.4).or.(j.eq.4.and.x(j).eq.3)) then
			!11 22 33 44 12 21 34 43
				kinetic3=2.d0**0.5d0*pi**1.5d0*miul*miur/(miul+miur)**2.5d0*(0.d0-0.d0)*Exp(-2.d0*miul*miur/(miul+miur)*(0.d0-0.d0))
			else
			!13 23 14 24 31 32 41 42
				kinetic3=2.d0**0.5d0*pi**1.5d0*miul*miur/(miul+miur)**2.5d0*(sl-sr)*Exp(-2.d0*miul*miur/(miul+miur)*(sl-sr))
			endif
			
			
			kinetic4=kinetic4+2.d0*kinetic2*kinetic3*overlap
		enddo
	enddo
kinetic=kinetic1-kinetic4
end


subroutine chazhi(y)
    implicit none
	real*8 xa(30),ya(30),x(1000),y(1000),xx,yy
    integer n,i
	open(2,file='30.30.dat',status='old')
		read(2,*) (xa(i),ya(i),i=1,30)
	close(2)
	n=30
	do i=1,1000
		xx=0.003*float(i)
		call polint(xa,ya,n,xx,yy)
		x(i)=xx
		y(i)=yy
	end do
	open(3,file='1000.1000.dat')
		do i=1,1000
			write(3,*) x(i),y(i)
		enddo
	close(3)
	return
end

SUBROUTINE POLINT(XA,YA,N,X,Y)
    implicit none
    integer i,n
    real*8 XA(30),YA(30),x,y
		do i=1,n
			if(abs(x-xa(i)).le.1.e-5) then
			y=ya(i)
			return
			endif
		enddo
		if(x.lt.xa(1)) then
			y=ya(2)+(x-xa(2))*(ya(1)-ya(2))/(xa(1)-xa(2))
		endif
		!if(x.gt.xa(30)) then
		!	y=0.d0
		!	goto end
		!endif
		DO I=1,29 
			if(x.gt.xa(i).and.x.lt.xa(i+1)) then
				y=ya(i)+(x-xa(i))*(ya(i+1)-ya(i))/(xa(i+1)-xa(i))
			endif
		enddo
      RETURN
 END
