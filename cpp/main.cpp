/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <string>
#include <iostream>
#include "Domain.h"
#include "Bar.h"
#include "FourQ.h"
#include "Outputter.h"
#include "Clock.h"
#include "GlobalsDefine.h"


using namespace std;

int main(int argc, char *argv[])//argc是命令行总的参数个数,argv[]是argc个字符串，其中第0个参数是程序的全名
{   
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
    }
    
	string filename(argv[1]);//用指针argv[1]所指向的字符串常量初始化string对象
    string InFile = filename + ".dat";
//  cout <<InFile << endl;
    string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::Instance();//调用静态函数，创建求解域类，并定义FEMData是求解域指针

    Clock timer;   
    timer.Start();//开启时钟

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))//读数据，并建立COutputter* Output
	
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();//返回时钟工作后所开启的时间（秒数）

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
//为力向量，列高，对角元地址和刚度矩阵分配空间
//计算列高度和对角线单元的地址 

	FEMData->AllocateMatrices();
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix(); 
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements//求解平衡方程
    CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
//对K矩阵分解以及用回代法求解位移
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

    COutputter* Output = COutputter::Instance();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)//对第LoadCase载荷类组装力矢量
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute//用回代法求解位移
        Solver->BackSubstitution(FEMData->GetForce());
            
#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
#endif            
        Output->OutputNodalDisplacement(lcase);
    }

    double time_solution = timer.ElapsedTime();
//  Calculate and output stresses of all elements
	Output->OutputElementStress();//计算应力并输出应力
    
    double time_stress = timer.ElapsedTime();
    
    timer.Stop();//关闭时钟
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	return 0;
}