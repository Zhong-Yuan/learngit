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

#include <string>
#include <fstream>
#include <vector>

#include "Node.h"
#include "Bar.h"
#include "ElementGroup.h"
#include "Outputter.h"
#include "Solver.h"
#include "LoadCaseData.h"
#include "SkylineMatrix.h"

using namespace std;

//!	Clear an array
template <class type> void clear( type* a, unsigned int N );

//!	Domain class : Define the problem domain
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:

//!	The instance of the Domain class
	static CDomain* _instance;//将变量存储在程序的静态存储区而非栈上空间，是所有对象所公有的

//!	Input file stream for reading data from input data file //ifstream默认是从硬盘到内存
	ifstream Input;

//!	Heading information for use in labeling the outpu
	char Title[256]; 

//!	Solution MODEX 求解模式
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	Total number of nodal points 节点总数
	unsigned int NUMNP;

//!	List of all nodes in the domain 域中节点列表 CNode类数组
	CNode* NodeList;

//!	Total number of element groups. 单元组数量
/*! An element group consists of a convenient collection of elements with same type */
//单元组是相同种类单元的集合（如杆单元集合等）
	unsigned int NUMEG;

	//! Element group list
    CElementGroup* EleGrpList;

////!	Element type of each group //单元组对应的单元种类（指针）
//	unsigned int* ElementTypes;
//
////!	Number of elements in each element group //单元组的单元数（指针）
//	unsigned int* NUME;
//
////!	Element Set List
///*!		ElementSetList[i] - ith element set */
///*!		ElementSetList[i][j] - jth element in ith set */
////单元设置列表 每行表示一个单元组设置
//	CElement** ElementSetList;
//
////!	Number of different sets of material/section properties in each element group 各单元组中的材料数目
//	unsigned int* NUMMAT;
//
////!	Material set list
///*!		MaterialSetList[i] - ith material set */
///*!		MaterialSetList[i][j] - jth material in ith set */
////材料设置列表 每行表示一个单元组材料设置
//	CMaterial** MaterialSetList;

//!	Number of load cases 载荷类数
	unsigned int NLCASE;

//!	List of all load cases 载荷类列表
	CLoadCaseData* LoadCases;

//!	Number of concentrated loads applied in each load case
//集中荷载的个数
	unsigned int* NLOAD; 

//!	Total number of equations in the system
	unsigned int NEQ;

//!	Number of elements in banded global stiffness matrix
//带状整体刚度矩阵中的单元个数
	unsigned int NWK;

//!	Maximum half bandwith
//最大半带宽
	unsigned int MK;

//!	Banded stiffness matrix
/*! A one-dimensional array storing only the elements below the	skyline of the 
	global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;
	//带状整体刚度矩阵的向量式

//!	Global nodal force/displacement vector
//全局力/位移矩阵
	double* Force;

public:

//!	Constructor
	CDomain();

//!	Desconstructor
	~CDomain();

//!	Return pointer to the instance of the Domain class
	static CDomain* Instance();

//!	Read domain data from the input data file
	bool ReadData(string FileName, string OutFile);

//!	Read nodal point data
	bool ReadNodalPoints();

//!	Read load case data
	bool ReadLoadCases();

//!	Read element data
	bool ReadElements();

//!	Read bar element data from the input data file
	bool ReadBarElementData(unsigned int EleGrp);

//!	Calculate global equation numbers corresponding to every degree of freedom of each node
	void CalculateEquationNumber();

//!	Calculate column heights
	void CalculateColumnHeights();

	void CalculateColumnHeightForCBar();

	void CalculateColumnHeightForFourQ();

//!	Calculate address of diagonal elements in banded matrix
	void CalculateDiagnoalAddress();

//! Allocate storage for matrices
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

//!	Assemble the banded gloabl stiffness matrix
	void AssembleStiffnessMatrix();

	void AssembleStiffnessMatrixForCBar();

	void AssembleStiffnessMatrixForFourQ();

//!	Assemble the global nodal force vector for load case LoadCase
	bool AssembleForce(unsigned int LoadCase); 
//为所有单元进行重力分配
	bool DisGForceall(double g,int Dg);

	bool BarDisGForce(unsigned int EleGrp,double g,int Dg);


//!	Return solution mode
	inline unsigned int GetMODEX() { return MODEX; }

//!	Return the title of problem
	inline string GetTitle() { return Title; }

//!	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//!	Return the total number of nodal points
	inline unsigned int GetNUMNP() { return NUMNP; }

//!	Return the number of banded global stiffness matrix elements
	inline unsigned int GetNWK() { return NWK; }

//!	Return the maximum half bandwith
	inline unsigned int GetMK() { return MK; }

//!	Return the node list
	inline CNode* GetNodeList() { return NodeList; }

////!	Return the number of elements in each element group
//	inline unsigned int* GetNUME() { return NUME; }
//
//!	Return total number of element groups
	inline unsigned int GetNUMEG() { return NUMEG; }

//! Return element group list
    CElementGroup* GetEleGrpList() { return EleGrpList; }

////!	Element type of each group
//	inline unsigned int* GetElementTypes() {return ElementTypes; }
//
////!	Return element Set List 
//	inline CElement** GetElementSetList() { return ElementSetList; }
//
////!	Return number of different sets of material/section properties in each element group
//	inline unsigned int* GetNUMMAT() { return NUMMAT; }
//
////!	Return material set list
//	inline CMaterial** GetMaterialSetList() { return MaterialSetList; }
//
//!	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//!	Return pointer to the global nodal displacement vector 这时力向量已被求解替换为位移向量
	inline double* GetDisplacement() { return Force; }

//!	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//!	Return the number of concentrated loads applied in each load case
	inline unsigned int* GetNLOAD() { return NLOAD; }

//!	Return the list of load cases
	inline CLoadCaseData* GetLoadCases() { return LoadCases; }

//!	Return pointer to the banded stiffness matrix
	inline CSkylineMatrix<double>* GetStiffnessMatrix() { return StiffnessMatrix; }

};
