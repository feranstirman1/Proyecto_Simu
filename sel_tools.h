node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        case 4: n = m.getNode(e.getNode4()-1); break;
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE:    v = n.getY(); break;
        case ZETA:  v = n.getZ(); break;
	}
	return v;
}

float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m), n2=selectNode(j,e,m);
	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float calculateMagnitude(float v1, float v2, float v3){
    return sqrt(pow(v1,2)+pow(v2,2)+pow(v3, 2));
}

float OperarRestaTenedor(element e, int coord1, int coord2,  float value_a, float value_b, mesh &m){
    float a, b, c, d;

    a = calcularTenedor(e, coord1, value_a, 1 ,m);
    b = calcularTenedor(e, coord2, value_b, 1 ,m);
    c = calcularTenedor(e, coord1, value_b, 1 ,m);
    d = calcularTenedor(e, coord2, value_a, 1 ,m);

    return (a*b)-(c*d);
}


float calcularLambda1(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = 3*(z1*z1) + 2*z1 + (z2+z3+z4) + (z2*z2) + z2*(z3+z4) + (z3*z3) + z3*(z4) + 18*(x1) + 9*(x2) + 9*(x3) + 9*(x4) + 12*(y1) + 6*(y2) + 6*(y3) + 6*(y4) + (z4*z4);

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularLambda2(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(2*z2 + z3 + z4 ) + 3*pow(z2,2) + 2*z2*(z3+z4) + pow(z3,2) + z4*(z3) + 9*(x1) + 18*(x2) + 9*(x3) + 9*(x4) + 6*y1 + 12*y2 + 6*y3 + 6*y4 + pow(z4,2);

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularLambda3(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(z2 + 2*z3 + z4 ) + pow(z2,2) + z2*(2*z3+z4) + 3*pow(z3,2) + 2*z4*(z3) + 9*(x1) + 9*(x2) + 18*(x3) + 9*(x4) + 6*y1 + 6*y2 + 12*y3 + 6*y4 + pow(z4,2);

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularLambda4(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(z2 + z3 + 2*z4 ) + pow(z2,2) + z2*(z3+z4) + pow(z3,2) + 2*z4*(z3) + 9*(x1) + 9*(x2) + 9*(x3) + 18*(x4) + 6*y1 + 12*y2 + 6*y3 + 12*y4 + 3*pow(z4,2);

	resultado = resultado / 360.0 ;

	return resultado;
}

float calcularGPrima1(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = 3*pow(z1,2) + 2*z1*(z2 + z3 + z4 ) + pow(z2,2) + z2*(z3+z4) + pow(z3,2) + z4*(z3) + 12*pow(x1,2) + 8*x1*(x2+x3+x4) + 4*pow(x2,2) + 4*x2*(x3+x4) + 4*pow(x3,2) + 4*x3*x4 + 4*pow(x4,2) + pow(z4,2); 

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularGPrima2(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(2*z2 + z3 + z4 ) + 3*pow(z2,2) + 2*z2*(z3+z4) + pow(z3,2) + z4*(z3) + 4*pow(x1,2) + 4*x1*(2*x2+x3+x4) + 12*pow(x2,2) + 8*x2*(x3+x4) + 4*pow(x3,2) + 4*x3*x4 + 4*pow(x4,2) + pow(z4,2); 

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularGPrima3(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(z2 + 2*z3 + z4 ) + pow(z2,2) + z2*(2*z3+z4) + 3*pow(z3,2) + 2*z4*(z3) + 4*pow(x1,2) + 4*x1*(x2+2*x3+x4) + 4*pow(x2,2) + 4*x2*(2*x3+x4) + 12*pow(x3,2) + 8*x3*x4 + 4*pow(x4,2) + pow(z4,2); 

	resultado = resultado / 360.0 ;

	return resultado;
}
float calcularGPrima4(element e,mesh &m){
	//se declara el resultado
	float resultado;
	
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	resultado = pow(z1,2) + z1*(z2 + z3 + 2*z4 ) + pow(z2,2) + z2*(z3+2*z4) + pow(z3,2) + 2*z4*(z3) + 4*pow(x1,2) + 4*x1*(x2+x3+2*x4) + 4*pow(x2,2) + 4*x2*(x3+2*x4) + 4*pow(x3,2) + 8*x3*x4 + 12*pow(x4,2) + 3*pow(z4,2); 

	resultado = resultado / 360.0 ;

	return resultado;
}


//funcion auxiliar para la matriz H
float calcularTermino(float var1,float var2,float var3,float var4, int termino){
	float resultado;
	switch(termino) {
  		case 1:
    		resultado = 3*pow(var1,2) + 2*var1*(var2 + var3 + var4 ) + pow(var2,2) + var2*(var3+var4) + pow(var3,2) + var4*(var3) + pow(var4,2);
    		break;
  		case 2:
    		resultado = pow(var1,2) + var1*(2*var2 + var3 + var4 ) + 3*pow(var2,2) + 2*var2*(var3+var4) + pow(var3,2) + var4*(var3) + pow(var4,2);
    		break;
		case 3:
    		resultado = pow(var1,2) + var1*(var2 + 2*var3 + var4 ) + pow(var2,2) + var2*(2*var3+var4) + 3*pow(var3,2) + 2*var4*(var3) + pow(var4,2);
    		break;
		case 4:
    		resultado = pow(var1,2) + var1*(var2 + var3 + 2*var4 ) + pow(var2,2) + var2*(var3+2*var4) + pow(var3,2) + 2*var4*(var3) + 3*pow(var4,2);
    		break;
  		default:
    		// code block
			break;
	}
	return resultado;
}

float calcularHPrima1(element e,mesh &m){
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	//se declaran los cuatro termino
	float a = calcularTermino(x1,x2,x3,x4,1) + calcularTermino(y1,y2,y3,y4,1) + calcularTermino(z1,z2,z3,z4,1);
	float b = calcularTermino(x1,x2,x3,x4,2) + calcularTermino(y1,y2,y3,y4,2) + calcularTermino(z1,z2,z3,z4,2);
	float c = calcularTermino(x1,x2,x3,x4,3) + calcularTermino(y1,y2,y3,y4,3) + calcularTermino(z1,z2,z3,z4,3);
	float d = calcularTermino(x1,x2,x3,x4,4) + calcularTermino(y1,y2,y3,y4,4) + calcularTermino(z1,z2,z3,z4,4);

	a = a / 360.0 ;
	b = b / 360.0 ;
	c = c / 360.0 ;
	d = d / 360.0 ;

	return a;
}
float calcularHPrima2(element e,mesh &m){
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	//se declaran los cuatro termino
	float a = calcularTermino(x1,x2,x3,x4,1) + calcularTermino(y1,y2,y3,y4,1) + calcularTermino(z1,z2,z3,z4,1);
	float b = calcularTermino(x1,x2,x3,x4,2) + calcularTermino(y1,y2,y3,y4,2) + calcularTermino(z1,z2,z3,z4,2);
	float c = calcularTermino(x1,x2,x3,x4,3) + calcularTermino(y1,y2,y3,y4,3) + calcularTermino(z1,z2,z3,z4,3);
	float d = calcularTermino(x1,x2,x3,x4,4) + calcularTermino(y1,y2,y3,y4,4) + calcularTermino(z1,z2,z3,z4,4);

	a = a / 360.0 ;
	b = b / 360.0 ;
	c = c / 360.0 ;
	d = d / 360.0 ;

	return b;
}
float calcularHPrima3(element e,mesh &m){
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	//se declaran los cuatro termino
	float a = calcularTermino(x1,x2,x3,x4,1) + calcularTermino(y1,y2,y3,y4,1) + calcularTermino(z1,z2,z3,z4,1);
	float b = calcularTermino(x1,x2,x3,x4,2) + calcularTermino(y1,y2,y3,y4,2) + calcularTermino(z1,z2,z3,z4,2);
	float c = calcularTermino(x1,x2,x3,x4,3) + calcularTermino(y1,y2,y3,y4,3) + calcularTermino(z1,z2,z3,z4,3);
	float d = calcularTermino(x1,x2,x3,x4,4) + calcularTermino(y1,y2,y3,y4,4) + calcularTermino(z1,z2,z3,z4,4);

	a = a / 360.0 ;
	b = b / 360.0 ;
	c = c / 360.0 ;
	d = d / 360.0 ;

	return c;
}
float calcularHPrima4(element e,mesh &m){
	//Se seleccionan los nodos a utilizar
	node n1=selectNode(1,e,m);
	node n2=selectNode(2,e,m);
	node n3=selectNode(3,e,m);
	node n4=selectNode(4,e,m);

	//se obtienen todos los valores que se van a utilizar
	float x1 = selectCoord(EQUIS,n1);
	float x2 = selectCoord(EQUIS,n2);
	float x3 = selectCoord(EQUIS,n3);
	float x4 = selectCoord(EQUIS,n4);

	float y1 = selectCoord(YE,n1);
	float y2 = selectCoord(YE,n2);
	float y3 = selectCoord(YE,n3);
	float y4 = selectCoord(YE,n4);

	float z1 = selectCoord(ZETA,n1);
	float z2 = selectCoord(ZETA,n2);
	float z3 = selectCoord(ZETA,n3);
	float z4 = selectCoord(ZETA,n4);

	//se declaran los cuatro termino
	float a = calcularTermino(x1,x2,x3,x4,1) + calcularTermino(y1,y2,y3,y4,1) + calcularTermino(z1,z2,z3,z4,1);
	float b = calcularTermino(x1,x2,x3,x4,2) + calcularTermino(y1,y2,y3,y4,2) + calcularTermino(z1,z2,z3,z4,2);
	float c = calcularTermino(x1,x2,x3,x4,3) + calcularTermino(y1,y2,y3,y4,3) + calcularTermino(z1,z2,z3,z4,3);
	float d = calcularTermino(x1,x2,x3,x4,4) + calcularTermino(y1,y2,y3,y4,4) + calcularTermino(z1,z2,z3,z4,4);

	a = a / 360.0 ;
	b = b / 360.0 ;
	c = c / 360.0 ;
	d = d / 360.0 ;

	return d;
}

float calculateLocalVolume(int i,mesh m){

    double Ve, u, v, w, U, V, W, a, b, c, d, X, x, Y, y, Z, z;
    
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    U = calculateMagnitude(n2.getX()-n1.getX(), n2.getY()-n1.getY(), n2.getZ()-n1.getZ() );
    V = calculateMagnitude(n3.getX()-n2.getX(), n3.getY()-n2.getY(), n3.getZ()-n2.getZ() );
    W = calculateMagnitude(n3.getX()-n1.getX(), n3.getY()-n1.getY(), n3.getZ()-n1.getZ() );

    u = calculateMagnitude(n4.getX()-n3.getX(), n4.getY()-n3.getY(), n4.getZ()-n3.getZ() );
    v = calculateMagnitude(n4.getX()-n1.getX(), n4.getY()-n1.getY(), n4.getZ()-n1.getZ() );
    w = calculateMagnitude(n4.getX()-n2.getX(), n4.getY()-n2.getY(), n4.getZ()-n2.getZ() );

    X = (w-U+v)*(U+v+w);
    x = (U-v+w)*(v-w+U);

    Y = (u-V+w)*(V+w+u);
    y = (V-w+u)*(w-u+V);

    Z = (v-W+u)*(W+u+v);
    z = (W-u+v)*(u-v+W);

    a = sqrt(x*Y*Z);
    b = sqrt(y*Z*X);
    c = sqrt(z*X*Y);
    d = sqrt(x*y*z);

    Ve = sqrt( (-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d) ) / (192*u*v*w);
    
    return Ve;

}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}
