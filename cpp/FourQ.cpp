

//平面四边形单元
#include "FourQ.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigen>
#include "GlobalsDefine.h"

using namespace std;

//	Constructor
FourQ::FourQ()
{
	NEN = 4;	// Each element has 4 nodes
	nodes = new CNode*[NEN];
    
    ND = 8;     //8个自由度
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
	Area = 0;    // 0面积是非法值，调用时可以检测
}

//	Desconstructor
FourQ::~FourQ()
{
	//不需要显式地调用基类的析构函数，系统会自动隐式调用。
	//先执行派生类析构函数的函数体，再调用基类的析构函数。
}

//	Read element data from stream Input 读取并配置单元组中第Ele号4Q单元的材料配置号（在对应材料组中的编号），左右节点号
bool FourQ::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
//MaterialSets实参是FourQMaterial* MaterialSetList[EleGrp]
{
	unsigned int N;
	Input >> N;	// element number
	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;
		//cerr：输出到标准错误的ostream对象，常用于程序错误信息
		return false;
	}
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial = &(dynamic_cast<FourQMaterial*>(MaterialSets))[MSet - 1];//FourQMaterial* ElementMaterial
	//dynamic_cast 将一个基类对象指针（或引用）cast到继承类指针
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];

	return true;
}
//	Write element data to stream
void FourQ::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) 
	       << nodes[0]->NodeNumber << setw(9) << nodes[1]->NodeNumber << setw(9) 
		   << nodes[2]->NodeNumber << setw(9) << nodes[3]->NodeNumber 
		   << setw(12) << ElementMaterial->nset << endl;
}

//!获得单元的面积
double FourQ::GetArea() 
{	
	//需要修改！
	return Area;
}
//!计算单元的面积
double FourQ::CalcArea() 
{
	//nodes[0]->XYZ[0],nodes[0]->XYZ[1], 1 3
	//nodes[1]->XYZ[0],nodes[1]->XYZ[1], 2
	//nodes[2]->XYZ[0],nodes[2]->XYZ[1], 3 1
	//nodes[3]->XYZ[0],nodes[3]->XYZ[1];   2
	//S=(1/2)*(x1y2+x2y3+x3y1-x1y3-x2y1-x3y2)
	//需要修改!考虑对此处进行代码优化
	Area = 0.5*(nodes[0]->XYZ[0]*nodes[1]->XYZ[1]+nodes[1]->XYZ[0]*nodes[2]->XYZ[1]+nodes[2]->XYZ[0]*nodes[0]->XYZ[1]
		   -nodes[0]->XYZ[0]*nodes[2]->XYZ[1]-nodes[1]->XYZ[0]*nodes[0]->XYZ[1]-nodes[2]->XYZ[0]*nodes[1]->XYZ[1]+
		   		nodes[2]->XYZ[0]*nodes[3]->XYZ[1]+nodes[3]->XYZ[0]*nodes[0]->XYZ[1]+nodes[0]->XYZ[0]*nodes[2]->XYZ[1]
		   -nodes[2]->XYZ[0]*nodes[0]->XYZ[1]-nodes[3]->XYZ[0]*nodes[2]->XYZ[1]-nodes[0]->XYZ[0]*nodes[3]->XYZ[1]);


	return Area;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
//  建立LM: 局部自由度对应的全局自由度 
void FourQ::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++) //每个单元节点数循环
	{
		for (unsigned int D = 0; D < 2; D++) //每个节点自由度循环(对4Q单元只使用节点的x和y值)
        {//注意到对于4Q LocationMatrix只有8个元素
			LocationMatrix[i++] = nodes[N]->bcode[D];//这时bcode存储的是全局自由度编号
			//cout << "Loca" << i << "=" << LocationMatrix[i-1] <<endl;
		}
	}

}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4Q element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int FourQ::SizeOfStiffnessMatrix() 
	{ 
		#ifdef _EIGEN_ //使用Eigen
			return 8; //矩阵为8x8
		#else //不使用Eigen
			return 36; 
		#endif
	}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
//  计算4Q单元的刚度阵
#ifdef _EIGEN_ //使用Eigen
void FourQ::ElementStiffness(CSkylineMatrix<double>* StiffnessMatrix)
{	
	Eigen::Matrix<double, 8, 8> Matrix= Eigen::Matrix<double, 8, 8>::Zero();
	FourQMaterial* material = dynamic_cast<FourQMaterial*>(ElementMaterial);	// Pointer to material of the element
	//clear(Matrix, SizeOfStiffnessMatrix());//清空矩阵
	int ngp=2;
	double gp[] = {-0.57735027, 0.57735027}; //用于两点Gauss积分    
    double w[]  = {1,          1          };
	double detJ;
	Eigen::Matrix<double, 4, 2> C;
	C << nodes[0]->XYZ[0],nodes[0]->XYZ[1],
		 nodes[1]->XYZ[0],nodes[1]->XYZ[1],
		 nodes[2]->XYZ[0],nodes[2]->XYZ[1],
		 nodes[3]->XYZ[0],nodes[3]->XYZ[1];
	for (unsigned int i = 0; i < ngp; i++) //Gauss插值点的数目
	{
       for (unsigned int j = 0; j < ngp; j++) 
	   {
		//eta = gp[i]; psi = gp[j]; 
        	Eigen::Matrix<double, 2, 8> N = NmatElast2D(gp[i],gp[j]);       // shape functions matrix 
			Eigen::Matrix<double, 3, 8> B = BmatElast2D(gp[i],gp[j],&C,&detJ);     // derivative of the shape functions 
        	Matrix = Matrix + w[i]*w[j]*B.transpose()*(material->D)*B*detJ;   // element conductance matrix 
        	//be = N*b(:,e);                     // interpolate body forces using element shape functions 
        	//fe = fe + w(i)*w(j)*N.transpose()*be*detJ;    // element nodal force vector     
	   }
	}
	for (unsigned int j = 0; j < ND; j++) 
	{
		unsigned int Lj = LocationMatrix[j];	// Global equation number corresponding to jthDOF of the element//单元j+1自由度对应	全局号Lj
		if (!Lj) 
			continue;
	//	Address of diagonal element of column j in the one dimensional element stiffness matrix
		for (unsigned int i = 0; i <= j; i++)
		{
			unsigned int Li = LocationMatrix[i];	// Global equation number corresponding to ith DOF of the element
			if (!Li) 
				continue;
			(*StiffnessMatrix)(Li,Lj) += Matrix(i,j);
			//cout << "Li="<<Li <<setw(8) <<"Lj=" << Lj << endl;
		}
	}
	//cout << "Matrix" << Matrix <<endl;
	//cout << endl;
	//cout << endl;
	return;
}
#else //不使用Eigen
void FourQ::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());//清空矩阵

}
#endif

//计算形函数矩阵的值
Eigen::Matrix<double, 2, 8> FourQ::NmatElast2D(double eta,double psi)
{
	Eigen::Matrix<double, 2, 8> N;
	//N  =  [N1   0    N2    0    N3   0    N4   0;     
    //  	 0    N1   0     N2   0    N3    0   N4];
	N << 0.25*(1-psi)*(1-eta), 0, 0.25*(1+psi)*(1-eta), 0, 0.25*(1+psi)*(1+eta), 0, 0.25*(1-psi)*(1+eta), 0,
		 0, 0.25*(1-psi)*(1-eta), 0, 0.25*(1+psi)*(1-eta), 0, 0.25*(1+psi)*(1+eta), 0, 0.25*(1-psi)*(1+eta);

	//N(0,0) = 0.25*(1-psi)*(1-eta); N(1,1) = 0.25*(1-psi)*(1-eta);
	//N(2,0) = 0.25*(1+psi)*(1-eta); N(3,1) = 0.25*(1+psi)*(1-eta);
	//N(4,0) = 0.25*(1+psi)*(1+eta); N(5,1) = 0.25*(1+psi)*(1+eta);
	//N(6,0) = 0.25*(1-psi)*(1+eta); N(7,1) = 0.25*(1-psi)*(1+eta);
	return N;
}

//计算形函数导数矩阵B的值
Eigen::Matrix<double, 3, 8> FourQ::BmatElast2D(double eta,double psi,Eigen::Matrix<double, 4, 2>* C,double* detJ)
{	
	Eigen::Matrix<double, 2, 4> GN;//梯度矩阵
	GN << 0.25*(eta-1), 0.25*(1-eta),  0.25*(1+eta), 0.25*(-eta-1),
		  0.25*(psi-1), 0.25*(-psi-1), 0.25*(1+psi), 0.25*(1-psi);
	Eigen::Matrix<double, 2, 2> J = GN*(*C);
	*detJ = J(0,0)*J(1,1)-J(0,1)*J(1,0);
	Eigen::Matrix<double, 2, 4> BB = J.inverse()*GN; //第1行对x求导 第2行对y求导
	Eigen::Matrix<double, 3, 8> B;
	//B = [ B1x      0     B2x     0      B3x    0      B4x     0  ; 
    //        0     B1y     0     B2y      0     B3y     0      B4y;  
    //      B1y     B1x    B2y    B2x     B3y    B3x    B4y     B4x];
	B << BB(0,0),       0, BB(0,1),       0, BB(0,2),       0, BB(0,3),       0,
		       0, BB(1,0),       0, BB(1,1),       0, BB(1,2),       0, BB(1,3),
		 BB(1,0), BB(0,0), BB(1,1), BB(0,1), BB(1,2), BB(0,2), BB(1,3), BB(0,3);
	return B;
}


//	Calculate element stress //4Q单元的应力计算
void FourQ::ElementStress(Eigen::Matrix<double, 3, 4>* stress, double* Displacement,double* XX,double* YY)
{	
	Eigen::Matrix<double, 3, 1> stress_temp;
	FourQMaterial* material = dynamic_cast<FourQMaterial*>(ElementMaterial);	// Pointer to material of the element
	double gp[] = {-0.57735027, 0.57735027}; //用于两点Gauss积分
	int ind = 0; //计数参数
	for (unsigned int i = 0; i < 2; i++)
	{
       for (unsigned int j = 0; j < 2; j++) 
	   {	
		   	double detJ;
			Eigen::Matrix<double, 4, 2> C;
			C << nodes[0]->XYZ[0],nodes[0]->XYZ[1],
				  nodes[1]->XYZ[0],nodes[1]->XYZ[1],
				  nodes[2]->XYZ[0],nodes[2]->XYZ[1],
				  nodes[3]->XYZ[0],nodes[3]->XYZ[1];
			Eigen::Matrix<double, 8, 1> de;
			de << 0,0,0,0,0,0,0,0;
			for (unsigned int i = 0; i < 8; i++)
			{
				if (LocationMatrix[i])
				de(i,0) = Displacement[LocationMatrix[i]-1];
			}
			//de << Displacement[LocationMatrix[0]-1],Displacement[LocationMatrix[1]-1],
			//	  Displacement[LocationMatrix[2]-1],Displacement[LocationMatrix[3]-1],
			//	  Displacement[LocationMatrix[4]-1],Displacement[LocationMatrix[5]-1],
			//	  Displacement[LocationMatrix[6]-1],Displacement[LocationMatrix[7]-1];

        	Eigen::Matrix<double, 2, 8> N = NmatElast2D(gp[i],gp[j]);       // shape functions matrix 
			Eigen::Matrix<double, 3, 8> B = BmatElast2D(gp[i],gp[j],&C,&detJ);     // derivative of the shape functions 
			//[X,Y]= [N(1,1) N(1,3) N(1,5) N(1,7)]*C; 
			XX[ind] = N(1,1)*C(0,0)+N(1,3)*C(1,0)+N(1,5)*C(2,0)+N(1,7)*C(3,0);
			YY[ind] = N(1,1)*C(0,1)+N(1,3)*C(1,1)+N(1,5)*C(2,1)+N(1,7)*C(3,1);
			//strain(:,ind) = B*de; 
            //stress(:,ind) = D*B*de;       // compute the stresses 
			//!!考虑优化
			stress_temp = (material->D)*B*de; 
			(*stress)(0,ind) = stress_temp(0,0);
			(*stress)(1,ind) = stress_temp(1,0);
			(*stress)(2,ind) = stress_temp(2,0);
			ind++;
	   }
	}
}