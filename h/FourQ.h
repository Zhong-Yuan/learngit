

#pragma once

#include "Element.h"
#include <Eigen/Eigen>

using namespace std;

//! 4Q element class
class FourQ : public CElement
{
public:

//!	Constructor
	FourQ();

//!	Desconstructor
	~FourQ();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);
//! 获得四边形的面积
	virtual double GetArea();
//! 计算四边形的面积
	virtual double CalcArea();

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	//virtual void ElementStiffness(double* Matrix);
	virtual void ElementStiffness(CSkylineMatrix<double>* StiffnessMatrix);

//!	Calculate element stress
	virtual void ElementStress(Eigen::Matrix<double, 3, 4>* stress, double* Displacement,double* X,double* Y);

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

//计算形函数矩阵的值
	virtual Eigen::Matrix<double, 2, 8> NmatElast2D(double eta,double psi);
//计算形函数导数矩阵B的值
	virtual Eigen::Matrix<double, 3, 8> BmatElast2D(double eta,double psi,Eigen::Matrix<double, 4, 2>* C,double* detJ);

private:
    double Area;//杨氏模量和材料密度等均在材料设置中
};
