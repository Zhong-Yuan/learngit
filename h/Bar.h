/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Bar element class
class CBar : public CElement
{
public:

//!	Constructor
	CBar();

//!	Desconstructor
	~CBar();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);
//!获得杆的长度
	virtual double GetLength();
//!计算杆的长度
	virtual double CalcLength();

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
#ifdef _EIGEN_ //使用Eigen
	virtual void ElementStiffness(CSkylineMatrix<double>* StiffnessMatrix);
#else
	virtual void ElementStiffness(double* Matrix);
#endif

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

private:
	double Len;//杆长度
};
