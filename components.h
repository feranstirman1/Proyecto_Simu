float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

void calculateMatrizAPrima(int i,Matrix &R,mesh m){
	zeroes(R,12,3);
    element e = m.getElement(i);
    float a,b,c,d;
    a = calcularLambda1(e,m);
    b = calcularLambda2(e,m);
    c = calcularLambda3(e,m);
    d = calcularLambda4(e,m);
    R.at(0).at(0) = a; R.at(0).at(1) = 0; R.at(0).at(2) = 0;
    R.at(1).at(0) = b; R.at(1).at(1) = 0; R.at(1).at(2) = 0;
    R.at(2).at(0) = c; R.at(2).at(1) = 0; R.at(2).at(2) = 0;
    R.at(3).at(0) = d; R.at(3).at(1) = 0; R.at(3).at(2) = 0;
    R.at(4).at(0) = 0; R.at(4).at(1) = a; R.at(4).at(2) = 0;
    R.at(5).at(0) = 0; R.at(5).at(1) = b; R.at(5).at(2) = 0;
    R.at(6).at(0) = 0; R.at(6).at(1) = c; R.at(6).at(2) = 0;
    R.at(7).at(0) = 0; R.at(7).at(1) = d; R.at(7).at(2) = 0;
    R.at(8).at(0) = 0; R.at(8).at(1) = 0; R.at(8).at(2) = a;
    R.at(9).at(0) = 0; R.at(9).at(1) = 0; R.at(9).at(2) = b;
    R.at(10).at(0) = 0; R.at(10).at(1) = 0; R.at(10).at(2) = c;
    R.at(11).at(0) = 0; R.at(11).at(1) = 0; R.at(11).at(2) = d; 
	
}

void calculateMatrizAro(int i,Matrix &R,mesh m){
	zeroes(R,12,3);
    element e = m.getElement(i);
	float a,b,c,d;
	a = calcularGPrima1(e,m);
    b = calcularGPrima2(e,m);
    c = calcularGPrima3(e,m);
    d = calcularGPrima4(e,m);
    R.at(0).at(0) = a; R.at(0).at(1) = 0; R.at(0).at(2) = 0;
    R.at(1).at(0) = b; R.at(1).at(1) = 0; R.at(1).at(2) = 0;
    R.at(2).at(0) = c; R.at(2).at(1) = 0; R.at(2).at(2) = 0;
    R.at(3).at(0) = d; R.at(3).at(1) = 0; R.at(3).at(2) = 0;
    R.at(4).at(0) = 0; R.at(4).at(1) = a; R.at(4).at(2) = 0;
    R.at(5).at(0) = 0; R.at(5).at(1) = b; R.at(5).at(2) = 0;
    R.at(6).at(0) = 0; R.at(6).at(1) = c; R.at(6).at(2) = 0;
    R.at(7).at(0) = 0; R.at(7).at(1) = d; R.at(7).at(2) = 0;
    R.at(8).at(0) = 0; R.at(8).at(1) = 0; R.at(8).at(2) = a;
    R.at(9).at(0) = 0; R.at(9).at(1) = 0; R.at(9).at(2) = b;
    R.at(10).at(0) = 0; R.at(10).at(1) = 0; R.at(10).at(2) = c;
    R.at(11).at(0) = 0; R.at(11).at(1) = 0; R.at(11).at(2) = d; 
}

float calculateL(int i,mesh m){
    element e = m.getElement(i);
	
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

	resultado = 2*pow(z1,2) + 2*z1*(z2 + z3 + z4 ) + 2*pow(z2,2) + 2*z2*(z3+z4) + 2*pow(z3,2) + 2*z4*(z3) + 5*(x1) + 5*(x2) + 5*(x3) + 5*(x4) + 5*y1 + 5*y2 + 5*y3 + 5*y4 + 2*pow(z4,2);

	resultado = resultado / 20.0 ;

	return resultado;
	
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD;
    float J,Determinant;
    
    /* [ A+K  G ]
       [  D   0 ]
    */

    //////////////////////////////////////////// Matrix A /////////////////////////////////////////////////////////
    /*
    	1/D * J * A' * alpha * Beta
    */
    Matrix matrizAPrima, Alpha, Beta;

    //calculando el determinante y el jacobiano
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);
	
    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(Determinant); //coeficiente de la matriz A
    calculateMatrizAPrima(e,matrizAPrima,m);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    // real_a(J/D) * A' * alpha * Beta 
    productRealMatrix(real_a, productMatrixMatrix(matrizAPrima,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);
    
    
    ///////////////////////////////////////////// Matrix K ////////////////////////////////////////////////////////////////////////
    /*
        (J * l) / D * Beta t * alpha t * alpha * beta 
    */
    Matrix Alpha_t,Beta_t;
	
    float L = calculateL(e,m);
	// l / D*D    
    float real_k = (float) (L)/(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);
	// (l / D*D) * Betat*alphaT * alpha *beta  
    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);

    
    ////////////////////////////////////// Matrix G //////////////////////////////////////////////////////////////////////////////////
    Matrix Omega,matrizAro;
	
	calculateMatrizAro(e,matrizAro,m);    
    float real_g = (float) (J/Determinant);

    calculateOmega(Omega);
    
    // (J/D) * matrizAro * alpha * omega
    productRealMatrix(real_g,productMatrixMatrix(matrizAro,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixG);

    //////////////////////////////////////////// Matrix D ////////////////////////////////////////////////////////////////////////////
    Matrix g_matrix,g_matrix_t,Omega_t;

    float real_d = (float)(J/(24*Determinant));
	
    calculateGamma(g_matrix);
    transpose(Omega, Omega_t);
    transpose(g_matrix,g_matrix_t);
    // (J/24*D) * omegaT * alphaT * gammaT
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);

    ///////////////////////////////////// Matrix M /////////////////////////////////////////////////////////////////////////////////////
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

void calculateVectorH(int i,Vector &h,mesh &m){
    element e = m.getElement(i);
    zeroes(h,4);
    h.at(0) = calcularHPrima1(e,m);
    h.at(1) = calcularHPrima2(e,m);
    h.at(2) = calcularHPrima3(e,m);
    h.at(3) = calcularHPrima4(e,m);
}

Vector createLocalb(int e,mesh &m){
    float J;
    element ele = m.getElement(e);
    
    Vector b,b_aux,f;
    Matrix g_matrix;
	
    calculateF(f, m);
    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    
    zeroes(b_aux,16);
    //Se calcula b pero con los ultimos 4 espacios siendo 0
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    // calular J*h
    
	b.at(12) = J*calcularHPrima1(ele,m);
    b.at(13) = J*calcularHPrima2(ele,m);
    b.at(14) = J*calcularHPrima3(ele,m);
    b.at(15) = J*calcularHPrima4(ele,m);
    
    return b;
}
