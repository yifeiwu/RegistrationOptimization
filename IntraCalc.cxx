#include <iostream>
#include <fstream>
#include <vector>
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkMath.h"
#include "vtkPointLocator.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkPointData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


#include <ctype.h>
#include <time.h>
#include "kdtree.h"



//Revision Notes
//
//-Added correspondences, now outputs node matches
//-Rotations and translations are now handled by matrix operations via Eigen. Outputs are stable for small translations and rotations. Tree is now pre-computed ahead of time.
//






using namespace Eigen;
using namespace std;


void set4x4_matrix(Matrix<double,3,3> a, vtkMatrix4x4 *matrix)
{
   int i, j;
   for (i=0; i<3; i++)
   {
	 for (j=0; j<3; j++)
	 {
	    matrix->SetElement(i,j,a(i,j));
	 }
	 
   }
}



void TransformPoints(double p[3], MatrixXd R, MatrixXd &myLRS_M){
int i;
for(i=0; i < myLRS_M.cols(); i++ ){

//cout<<"pretrans :"<<myLRS_M.col(i).transpose()<<endl;

	myLRS_M.col(i)=R*myLRS_M.col(i);
	myLRS_M(0,i)=myLRS_M(0,i)+p[0];
	myLRS_M(1,i)=myLRS_M(1,i)+p[1];
	myLRS_M(2,i)=myLRS_M(2,i)+p[2];
	
	
//cout<<"posttrans :"<<myLRS_M.col(i).transpose()<<endl;

	
	
	}
}//End TransformPoints


void Objective_Function (VectorXd &fvec,VectorXd x, MatrixXd myBelm_M, MatrixXd myLRS_M, kdtree *ptree)
{

//build Rotation and Translation matrix
 
double theta[3];
double p[3];

p={x[0],x[1],x[2]};
theta={x[3],x[4],x[5]};

MatrixXd R(3,3);

R <<  0,-theta[2],theta[1],
      theta[2],0,-theta[0],
     -theta[1],theta[0],0;
     

R= R.exp();

//cout<<R<<endl;
//cout<<"p: "<<p[0]<<p[1]<<p[2]<<endl;

 int i=0;
 int closest_node_id;

 double pt[3], pt2[3];
 double shift_magn=0.0;
 double shift_vec[3];
 shift_vec[0]=0.0; shift_vec[1]=0.0; shift_vec[2]=0.0;


//Perform the transform

TransformPoints(p, R, myLRS_M);
 
//Build Locator

//kdtree *ptree; //Allocate structures

int *pch;
struct kdres *presult;



//for(i=0; i < myBelm_M.cols(); i++ ){
//kd_insert3(ptree, myBelm_M(0,i), myBelm_M(1,i) , myBelm_M(2,i), 0);
//}



//Calc displacement distance

for(i=0; i < myLRS_M.cols(); i++ ){

pt = { myLRS_M(0,i), myLRS_M(1,i) , myLRS_M(2,i)};
presult = kd_nearest(ptree,pt);
pch = (int*)kd_res_item( presult, pt2 );

cout<<"Node Correspondence "<<i<<" -> "<<*pch<<endl;


   shift_vec[0] = pt[0]-pt2[0];
   shift_vec[1] = pt[1]-pt2[1];
   shift_vec[2] = pt[2]-pt2[2];
   
   
//   cout<<"shift"<<shift_vec[0]<<endl;
//   cout<<"shift2"<<shift_vec[1]<<endl;
//   cout<<"shift3"<<shift_vec[2]<<endl;
//   
   fvec(i*3)   =  shift_vec[0];
   fvec(i*3+1) =  shift_vec[1];
   fvec(i*3+2) =  shift_vec[2];

}
}//end Objective Function
    
    
    
void Calc_Jacobian(int num_nodes, int n_par,int m_dat, VectorXd params, MatrixXd &J, MatrixXd Belm_M, MatrixXd LRS_M, kdtree *ptree){

cout<<"calc jacobian"<<endl;
double delta = 1e-6;
VectorXd res(3*num_nodes);
VectorXd temp_res(3*num_nodes);
VectorXd temp_params;


Objective_Function(res, params, Belm_M, LRS_M, ptree);


int i;
for (i=0;i<n_par;i++)
        {
        temp_params = params;
        
        temp_params(i)=params(i)+delta;
//        cout<<temp_params.transpose()/delta<<endl;
        Objective_Function(temp_res,temp_params,Belm_M, LRS_M,ptree);
           
//        cout<<"res "<<res<<endl;
//        cout<<"temp res "<<temp_res<<endl;
        J.col(i)=(temp_res-res)/delta;
//        cout<<"jcol "<<i<<" "<<J.col(i)<<endl;

    }

}//end Calc_Jacobian



    
    




int main( int argc, char * argv[] ){


 int i, j, k, l;

 // get the base directory for reading and writing files 
 // obtaining it as a string makes life easier when writing files
 string base_directory = "/home/caleb/Desktop/janet/AllAnt_noPost_TRUTH_M/";

 // case id
 string caseID = "03282012";

 // get the patient id  
 string patientID;  

 patientID = "0328";  

 // get the pre-operative file names
 string preop_vtk_bel_filename = base_directory+"PreOperative/"+patientID+"_bel.vtk";
 string preop_nml_bel_filename = base_directory+"PreOperative/"+patientID+".bel";
 string preop_nml_nod_filename = base_directory+"PreOperative/"+patientID+".nod";

 // get the intra-op file names
// string intraop_xformed_srf_filename = base_directory+"IntraOperative/"+caseID+"_srf.vtk";
string intraop_xformed_srf_filename = base_directory+"IntraOperative/"+"phantomcleaned.vtk";

 // number of nodes in the volumetric mesh and number of bdry elements
 int number_of_nodes;
 int number_of_BdryElms;
 
 cout << "Reading Inputs" << endl;

 // read the bdry elm polydata with flags and the registered scanner surface: use these datasets to compute closest point distances
 vtkPolyDataReader * BdryElmReader = vtkPolyDataReader::New();
 BdryElmReader->SetFileName( preop_vtk_bel_filename.c_str() );
 BdryElmReader->Update();

 // read the registered scanner surface
 vtkPolyDataReader * LRSReader = vtkPolyDataReader::New();
 LRSReader->SetFileName( intraop_xformed_srf_filename.c_str() );
 LRSReader->Update();

 // get the polydatas for bdry elements and registered scanner surface
 vtkPolyData * BdryElms = vtkPolyData::New();
 BdryElms->DeepCopy( LRSReader->GetOutput() );

 vtkPolyData * LRS = vtkPolyData::New();
 LRS->DeepCopy( LRSReader->GetOutput() );
 
//Done Reading, making data structures



//Convert LRS from mm to m *1e-3, store as Eigen Matrix

double pt[3];
//int num_LRS_pts = LRS->GetNumberOfPoints();
//int num_Belm_pts = BdryElms->GetNumberOfPoints();

int num_LRS_pts = 155;
int num_Belm_pts = 155;

MatrixXd LRS_M(3,num_LRS_pts);
MatrixXd Belm_M(3,num_Belm_pts);


 // loop thru all closest points on lrs 
for ( i=0; i < num_LRS_pts; i++ )
   {
   // get the position
   LRS->GetPoint( i,pt );


   LRS_M(0,i) =  pt[0];
   LRS_M(1,i) =  pt[1];
   LRS_M(2,i) =  pt[2];

   }

for ( i=0; i < num_Belm_pts; i++ )
   {
   // get the position
   BdryElms->GetPoint( i,pt );

 
   
   Belm_M(0,i) =  pt[0]+5;
   Belm_M(1,i) =  pt[1]+.1;
   Belm_M(2,i) =  pt[2]+0;

   }

double theta[3];
double p[3];
//theta[0]=.003;
//theta[1]=0.001;

//Matrix<double,3,3> R_init;

//R_init <<  0,-theta[2],theta[1],
//      theta[2],0,-theta[0],
//     -theta[1],theta[0],0;

//R_init= R_init.exp();


//TransformPoints(p, R_init, LRS_M);

int data[num_Belm_pts];

//Create the KD Tree for the Belm set
kdtree *ptree; //Allocate structures
ptree =kd_create(3);
cout<<"hi"<<endl;
for(i=0; i < Belm_M.cols(); i++ ){
data[i] =  i;
 cout<<"insertion err"<<endl;
kd_insert3(ptree, Belm_M(0,i), Belm_M(1,i) , Belm_M(2,i), &data[i] );

cout<<&i<<endl;
cout<<i<<endl;

}



// // testing
//for ( i=0; i < 3; i++ )
//   {
//   // get the position


//   LRS_M(0,i) =  1;
//   LRS_M(1,i) =  0;
//   LRS_M(2,i) =  i;

//   }

//for ( i=0; i < 3; i++ )
//   {
//   // get the position

// 
//   Belm_M(0,i) =  1;
//   Belm_M(1,i) =  0;
//   Belm_M(2,i) =  i;

//   }








//cout<<setprecision(11)<<LRS_M.col(0)<<endl;
//cout<<setprecision(11)<<LRS_M(0,0)+1.5e-10<<endl;

 



  //compute the minimum
  int m_dat;
  int num_nodes = num_LRS_pts;
  m_dat=num_nodes*3;
  int n_par =6;
  
  int info;
  double fnorm, covfac;
  VectorXd x(n_par);
  VectorXd x_delta(n_par);
  /* the following starting values provide a rough fit. */
  x.setConstant(n_par, 0.);
  
//  x<< 
//   -0.0143,
//    0.0055,
//   -0.0036,
//    0.0647,
//    0.0742,
//   -0.0555;
 

  // do the computation
  
  MatrixXd J(m_dat,n_par);
  MatrixXd A(4,4);
  MatrixXd JTJ, d, D;
  
  
  

  VectorXd res(3*num_nodes);

    double lambda=0;
    int num_iter=20;
    VectorXd b;
    


    for (i=0;i<num_iter;i++)
    {
    

        Calc_Jacobian(num_nodes, n_par, m_dat, x, J, Belm_M, LRS_M, ptree);
        Objective_Function(res, x, Belm_M, LRS_M, ptree);


//        A=J.transpose()*J+lambda*(MatrixXd::Identity(n_par,n_par));
        
        JTJ=J.transpose()*J;
        d=JTJ.diagonal();
        D=d.asDiagonal();
        A=JTJ+lambda*D;
        
        
//        cout<<"J="<<J<<endl;
        b=-J.transpose()*res;
        x_delta=A.inverse()*b;
        x=x+x_delta;
//        cout<<"x"<<x<<"xdelta "<<x_delta<<endl;
        cout<<"solution"<<endl;
        std::cout << x.transpose() << std::endl;
    }



//Final Transform



p={x[0],x[1],x[2]};
theta={x[3],x[4],x[5]};


Matrix<double,3,3> R_final;

R_final <<  0,-theta[2],theta[1],
      theta[2],0,-theta[0],
     -theta[1],theta[0],0;

R_final= R_final.exp();

//Define transform 
vtkMatrix4x4 *tmatrix = vtkMatrix4x4::New();
set4x4_matrix(R_final, tmatrix);
tmatrix->SetElement(0,3,p[0]);
tmatrix->SetElement(1,3,p[1]);
tmatrix->SetElement(2,3,p[2]);
tmatrix->SetElement(3,3,1.0);

/////////



vtkTransform * fxform = vtkTransform::New();
fxform->SetMatrix( tmatrix );

 vtkTransformPolyDataFilter * fxformer = vtkTransformPolyDataFilter::New();
 fxformer->SetInput(LRS );
 fxformer->SetTransform( fxform );
 fxformer->Update();
 LRS->DeepCopy( fxformer->GetOutput() );
  
//  fxform->Scale(1e3,1e3,1e3);
//  fxformer->SetInput( LRS_post );
// fxformer->SetTransform( fxform );
// fxformer->Update();
//  vtkPolyData * LRS_final = vtkPolyData::New();
//  LRS_final->DeepCopy( fxformer->GetOutput() );
//  
  

  vtkPolyDataWriter * writer = vtkPolyDataWriter::New();

  writer->SetFileName( "test.vtk" );
  writer->SetInput( LRS);
  writer->Write(); 







//clear memory




 return 0; // success
}


