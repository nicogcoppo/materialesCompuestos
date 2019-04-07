
######### CONFIGURACION ##########

pkg load symbolic;

warning ('off','OctSymPy:sym:rationalapprox');

msj=0;save ./octaveLamIt msj;save ./octaveLamSo msj;

######### DECLARACIONES #######################

T=sym(zeros(3,3));

syms O;

######### CARACTERISTICAS MEC. MATERIAL (MKS) ##########

E1=134e6;

E2=7000e6;

V21=0.25;

G12=4200e6;

V12=(E1/E2)*V21;

sigmaMaxTraccion=[1270e6 42e6 63e6]; # X Y xy

sigmaMaxCompresion=[-1130e6 -141e6 -63e6];

########## CARACTERISTICAS DEL LAMINADO ################

VF=0.6;

AW=0.145;

n=1;

rho=1530;

t=(AW*n)/(rho*VF);

nCAPAS=2;

########## MATRIZ DE CARGAS Y DE MOMENTOS ############################

N=[1000 50 0];

M=[0 5 0];

########### MATRIZ CONSTITUTIVA ########################

S(1,1)=1/E1;

S(1,2)=-V21/E2;

S(2,1)=-V12/E1;

S(2,2)=1/E2;

S(3,3)=1/G12;

########### MATRIZ DE RIGIDEZ ##########################

Q=inv(S);

########## MATRIZ DE TRANSFORMACION ####################

T(1,1)=cos(O)^2;

T(2,2)=cos(O)^2;

T(1,2)=sin(O)^2;

T(2,1)=sin(O)^2;

T(3,1)=-cos(O)*sin(O);

T(2,3)=-2*cos(O)*sin(O);

T(1,3)=2*cos(O)*sin(O);

T(3,2)=cos(O)*sin(O);

T(3,3)=cos(O)^2-sin(O)^2;

########## MATRIZ DE CONSTANTES #########################

R(1,1)=1;

R(2,2)=1;

R(3,3)=2;

##-- AJUSTES DEL METODO ITERATIVO

DIVISIONES=3; # Subdivisiones del rango de alfa

RANGO=90; # Angulo maximo para alfa en grados

##--

alfa=linspace(0,RANGO*pi/180,DIVISIONES);

CONVERGENCIA=0;

##- VECTOR ITERACIONES 

IT=zeros(1,nCAPAS*DIVISIONES);
for i=1:nCAPAS
  for j=1:DIVISIONES
    IT(1,j+(i-1)*DIVISIONES)=j;
  endfor
endfor

##-- ITERACIONES

soluciones=zeros(1,5+nCAPAS);

while CONVERGENCIA<1

  MPO=unique(nchoosek(IT,nCAPAS),'rows'); # Tomo filas distintas unicamente.
  
  it=0;

  so=0;
  
  for i=1:size(MPO,1)

    TEST=MPO(i,:); # Tomo de a una forma de orientacion
    
    try

      QM=0;

      A=0;

      D=0;

      B=0;
      
      for j=1:size(TEST,2) # CALCULO DE LA MATRIZ A y D
	
	Tn=double(subs(T,O,alfa(TEST(j))));  

	Qi=inv(Tn)*Q*R*Tn*inv(R); 

	try

	  QM=[QM,Qi];  
	  
	catch

	  QM=Qi;

	end_try_catch  
	
	A=A+Qi; # MATRIZ DE RIGIDEZ PLANA DEL LAMINADO

	zj=j*t-t*nCAPAS/2;

	zj_1=(j-1)*t-t*nCAPAS/2;

	D=D+(1/3)*Qi*(zj^3-zj_1^3); # MATRIZ DE RIGIDEZ A LA FLEXION

	B=B+(1/2)*Qi*(zj^2-zj_1^2); # MATRIZ DE CONECTIVIDAD
	
	
      endfor

      A=A*t;

      IA=inv(A);

      #-- MATRICES [a] [b] [d]

      Bm=-IA*B;

      Cm=B*IA;

      Dm=D-(B*IA)*B;

      d=inv(Dm);
      
      a=IA-(Bm*d)*Cm;

      b=Bm*d;
 
      #-- CALCULO DE LA DEFORMACION POR TENSION PLANA

      e0=a*N'+b*M';

      #-- CALCULO DE LA DEFORMACION POR FLEXION

      K=b*N'+d*M';

      #keyboard;

      #-- CALCULO DE LAS TENSIONES EN CADA LAMINADO
     
      for j=1:size(TEST,2)

	z=-(nCAPAS^2+(1-2*j)*nCAPAS)*t/(2*nCAPAS-2);

	eTotal=(e0+z*K);
	
	on=QM(1:3,j*3-2:j*3)*eTotal;

		
	if(on<=sigmaMaxTraccion' && on>=sigmaMaxCompresion')

	  FLAG(j)=1;

	else

	  FLAG(j)=0;

	endif
	
	
	
      endfor

      defMax=(e0+nCAPAS*t*K/2); # DEFORMACION MAXIMA LAMINADO LITERATURA

      defMaxProm=sum(defMax)/3; # Promedio de la deformacion maxima en los 3 grados de libertad
      
      vistaAngulos=sprintf('%i ',alfa(TEST)*180/pi);

      per=((so+it+1)/size(MPO,1))*100;

      #keyboard
      
      if(sum(FLAG)==size(FLAG,2))
		
	so++;

	try

	  soluciones=[soluciones;nCAPAS,defMaxProm,alfa(TEST)*180/pi,defMax'];
	  
	catch

	  soluciones=[nCAPAS,defMaxProm,alfa(TEST)*180/pi,defMax'];

	end_try_catch
	
	mejorIntento=sortrows(soluciones,[2])(1,:);
	
	printf ("Iteraciones: %i \n",it);

	printf ("Soluciones: %i \n",so);

	printf ("Mejor solucion -> CAPAS=%d MAX.DEF.=%d ANGULOS=%s \n",mejorIntento(1),mejorIntento(2),vistaAngulos);

	printf ("--------%d %%-------\n\n",per);

	
	[fi,c]=size(soluciones);

	msj=soluciones(fi,:);

	save ./octaveLamSo -append msj 

	
      else
	
	it++;

	try

	  iteracion=[iteracion;nCAPAS,defMaxProm,alfa(TEST)*180/pi,defMax'];
	  
	catch

	  iteracion=[nCAPAS,defMaxProm,alfa(TEST)*180/pi,defMax'];

	end_try_catch
	
	mejorIntento=sortrows(iteracion,[2])(1,:);
	
	printf ("Iteraciones: %i \n",it);

	printf ("Soluciones: %i \n",so);
	
	printf ("Mejor solucion -> CAPAS=%d MAX.DEF.=%d ANGULOS=%s \n",mejorIntento(1),mejorIntento(2),vistaAngulos);
	
	printf ("--------%d %%-------\n\n",per);

	[fi,c]=size(iteracion);

	msj=iteracion(fi,:);

	save ./octaveLamIt -append msj 


      endif
        				   
     end_try_catch
           
  endfor

#% RECONFIGURACION DEL PROCESO ITERATIVO



if (size(soluciones)(1)<2) 
    
    
  # SI SE PRECISAN MAS CAPAS
  
  nCAPAS++;
  
  ##- VECTOR ITERACIONES 

  IT=zeros(1,nCAPAS*DIVISIONES);
  for i=1:nCAPAS
    for j=1:DIVISIONES
      IT(1,j+(i-1)*DIVISIONES)=j;
    endfor
  endfor

else ## SI SE PRECISA OPTIMIZAR UNA SOLUCION HALLADA

  #keyboard
  
  DIVISIONES++;
  
  alfa=linspace(0,RANGO*pi/180,DIVISIONES);

  IT=zeros(1,nCAPAS*DIVISIONES);
  for i=1:nCAPAS
    for j=1:DIVISIONES
      IT(1,j+(i-1)*DIVISIONES)=j;
    endfor
  endfor

  
#################### CHEQUEO DE CONVERGENCIA

  try

    CONVE=sortrows(soluciones,[2])(1:2,2);
    
    if (CONVE(2)/CONVE(1)<=1.05)

      CONVERGENCIA=1;

    else

      CONVERGENCIA=0;

    endif

  catch

    CONVERGENCIA=0;

  end_try_catch
  
    
endif

endwhile

mejorIntento=sortrows(soluciones,[2])(1,:);

vistaAngulos=sprintf('%i ',mejorIntento(1,3:2+nCAPAS));

vistaDeformaciones=sprintf('%i ',mejorIntento(1,3+nCAPAS:size(soluciones,2)));

printf ("\n\nCONVERGENCIA ALCANZADA \n",it);

printf ("SOLUCION -> CAPAS=%d MAX.DEF.=%d ANGULOS=%s DEFORMACIONES=%s \n",mejorIntento(1),mejorIntento(2),vistaAngulos,vistaDeformaciones);

printf ("--------%d %%-------\n\n",per);
