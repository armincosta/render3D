/**
*	Tomolog3DMatching && Visualization (3D-Shape1, CT, 3D-Shape2)
*
*	
*	Author: Armin Costa
*	(c)opyright 2006-2010
*	
*	e-mail: armincosta@hotmail.com
*
*
*	This module implements a 3D rendering algorithms to match heterogenuous Log data.
*	It renders DiShape 3D data of the logs with and without bark. Note that this
*	might not be an easy task as the logs in bark can have debarked parts.
*	Once the matching has occurred, the XRay projection is also rendered as volume data. In this way
*	We get a true correlation between the 3 different measures. Huge chunks of .mif files can be
*	analyzed and results be print to a file.
*
*	file: render3D.cpp
*
*	This source code is distributed under the GNU GENERAL PUBLIC LICENSE Version 3 in the hope it might be helpful
*
*	
**/
#ifndef __RENDER_HEADER__
#define __RENDER_HEADER__

#include "..\lib\SLXData\shapedata.h"
#include "..\lib\SLXData\XRayStructs.h"
#include "..\lib\SLXData\XRayData.h"
#include "..\lib\SLXData\PhotoData.h"
#include "..\Core\LogInfo3D.h"
#include "..\lib\SLXData\XRayNumeric.h"

//include the INTEL numeric libraries
#define             nsp_UsesAll         //to include all the possible nsp headers
#include "..\..\lib\nsp\nsp.h"



/*
#include <Magick++.h>				// we use Image Libs from ImageMagick to operate on image formats
using namespace Magick;
*/


#include <windows.h>
#include <gl\gl.h>
#include <gl\glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define VISUALIZE_TOMO_GEOMETRY 1
#define LOAD_VISUAL_XRAY_CALIB 1 // enable this to load and visualize raw xray .mif files, useful for Calibration
								// set also LOAD_ONLY_XRAY = 1 ,     uses getXRAYSliceInfoByNumber



//#####################################################
// Directory Settings for Data and Data-Analysis Files
#define DIR_HASSLACHER    1
#define DIR_LOCAL_LINZ    2


#define CURRENT_DIR_SETTINGS  DIR_LOCAL_LINZ 
//#####################################################


//############# DEFINES ######################################################################
#define COLLECT_DATA 0 // toogle this if you want to use pure Data Evaluation without 3D Graphics && MSG's
#define MAX_MIF_FILES 15000   // we want to analyze max 10000 logs    iff COLLECT_DATA == 1

#define KR_ANALYSIS 0 //  switch this if you want export Log Curvature Data instead of true 3D Log data  USE with  COLLECT_DATA == 1
// iff KR_ANALYSIS == 1  {
#define SIMULATE_WZLR_SIGNAL 0  //  switch this for Dataset where mir->reserve[0]  in not used (Wurzel-Reduzierer) 
#define EXPORT_WA_PHOTOS 1  // this exports the pics of the Logs having WA == 1  to LOG_PHOTO_DIR\\%mif_file.mif.jpg
#define USE_PERCENTAGE_DIAM_RATION 0    // this calculates the Krümmung as stated in "Österreichische Holzhaldelsusances 1973"
#define LOAD_DISHAPE3D_MODELS 0
#define SMOOTH_WIN_SIZE 15   // 20 cm smoothing window
#define ALWAYS_TAKE_3D1_MEASURE 0
//  }

#define LOAD_ONLY_XRAY 1  // this loads raw Tomolog .mif objects, to analyze "Prüfkörper" data
#define CALC_XRAX_BARK_DETECT 0	
#define	EXPORT_DENSITY_PROFILES 0   // if LOAD_ONLY_XRAY == 1

//##############################################################################################
// Data File definitions
#if CURRENT_DIR_SETTINGS == DIR_LOCAL_LINZ
	#define ANALYSIS_DATA_FILE  "C:\\log3D.csv"

	#define ANALYSIS_DATA_FILE_KR  "C:\\log3DKr.csv"
	#define LOG_PHOTO_DIR "C:\\LOG_PICS"

	#define ANALYSIS_DATA_FILE_EICH "C:\\log3DEich.csv"
	#define ANALYSIS_DATA_FILE_DENSITY "C:\\Documents and Settings\\cos\\Eigene Dateien\\WORKDOCS\\TOMO_Hasslacher\\TestsHFA\\log3DDensity.csv"

	#define MIF_FILEDIR  "E:\\SLX_DATA"  
	#define CT_FILES "E:\\CTData\\"

#elif CURRENT_DIR_SETTINGS == DIR_HASSLACHER

	#define ANALYSIS_DATA_FILE  "C:\\log3D.csv"

	#define ANALYSIS_DATA_FILE_KR  "C:\\log3DKr.csv" 
	#define LOG_PHOTO_DIR "C:\\LOG_PICS"

	#define ANALYSIS_DATA_FILE_EICH "C:\\log3DEich.csv"
	#define ANALYSIS_DATA_FILE_DENSITY "C:\\log3DDensity.csv"

	#define MIF_FILEDIR  "F:\\datenKruemmung"// "E:\\SLX_DATA"  // "F:\\MIF_FILES"   //"F:\\datenKruemmung"
	#define CT_FILES "C:\\CTData\\"

#endif
//##############################################################################################

//##############################################################################################
// ONLY for Krümmung  iff KR_ANALYSIS == 1 
#define MIN_EINSEITIGE_KR 30 //min Diff. einseitige Kr
#define OFFSET_WURZEL_ANLAUF 5  // we start 5 cm from the log beginning
#define OFFSET_WURZEL_ANLAUF_B (OFFSET_WURZEL_ANLAUF+25)  // b section
#define OFFSET_WURZEL_ANLAUF_C (OFFSET_WURZEL_ANLAUF_B+60)
#define OFFSET_WURZEL_ANLAUF_D (OFFSET_WURZEL_ANLAUF_B+25)
#define DIFF_OFFSET_WURZEL_ANLAUF 15

#define POS_AVG_OFFSET 30   // in mm
#define MAX_WZLR_DIFF 250   // diff in 1/10 mm
//##############################################################################################


//#############  Defines for XRAY 3D rendering ##########################################
#define CT_IMG_PIX_SIZE 768
#define ACIS_THR 100
#define GRID_CELL_SIZE 1  // CT_IMG_PIX_SIZE/GRID_CELL_SIZE Factor
#define BROKEN_PIX_THRES  600
//#######################################################################################



//#############	 Defines for XRAY Calibration
#define MIN_PIXEL_REST 4				// this is used to control Measure overrun
#define TOTAL_NR_PIXELS 768
//#######################################################################################
#define	LOG_RADIAL_STEP	3   // optimal 3

// this is the reference measure  degree (TOMOLOG)
#define DEGREE_MEASURE 39   // D1
#define DEGREE_MEASURE_D1 39   // D1
#define DEGREE_MEASURE_D2  (180-DEGREE_MEASURE)

#define LENGTH_DETAIL 2     // optimal 2
#define SAMPLE_INTERVAL 20

#define ROTATION_FAKTOR 2    // we turn in a 1 degree step when matching logs 

#define USE_OFFSET 0
#define TOMOLOG_WB    576   //mm
#define TOMOLOG_WOB   470   // mm

#define TOMOLOG_WB2    0   //mm     it seems that only Sensor1 has this offset
#define TOMOLOG_WOB2   0   // mm

#define bClambDiShapeDiam false

#if COLLECT_DATA
#define USE_GUI_FILESELECT 0
#else
#define USE_GUI_FILESELECT 1
#endif

#define USE_MULTIPLE_STEP_MATCHING 0

#define TRUST_SLX_MATCHING 1	// enable this if the logs have already been matched by SLX (of course we do!!:), in this case the mached log is at index+1

#if COLLECT_DATA
#define MAX_LOGS3D 1		// maximal nr. or simultaneous loadable objects
#else
#define MAX_LOGS3D 10
#endif


//#########################################################################################
#define THREAD_TIMEOUT  20000

#define SQR(x)  ((x)*(x))

enum TP  {DISHAPE_WBARK, XRAY, DISHAPE_WOBARK};


class Model3D;  // forward declaration

long WINAPI MatchFunc(Model3D *obj);


#define MSG(x) if(!COLLECT_DATA){   MessageBox(NULL, x, "Msg", MB_OK);   }  \



// this struct ist used to mark line segments in Kr. Analysis
typedef struct tagKR_LINES {
	float m_1_x;
	float m_1_y;
	float m_1_z;
	float m_2_x;
	float m_2_y;
	float m_2_z;

	float o_1_x;
	float o_1_y;
	float o_1_z;
	float o_2_x;
	float o_2_y;
	float o_2_z;

	bool krAvail;

}KR_LINES;


  // Helper functions
static double radMath2GradUhr( double w );
static double gradUhr2RadMath(double w);
static double radMath2Grad( double w );
static double Grad2radMath(double w);
static double Grad2Uhr( double w );
static double Uhr2Grad( double w );



/*
	This class && Interface enables threading functionality for a given 3D Model
  */
struct IRunnable {
  virtual void run() = 0;
};


class Model3DThread {
public:
	Model3DThread(){

	}
  Model3DThread(IRunnable *ptr) {
    _threadObj = ptr;
  }
  ~Model3DThread(){
	 stopThread__();
  }
  void startThread__() {
    DWORD threadID;
    hThread = ::CreateThread(0, 0, threadProc, _threadObj, 0, &threadID);
	SetThreadPriority (hThread, THREAD_PRIORITY_NORMAL);   //THREAD_PRIORITY_LOWEST   THREAD_PRIORITY_NORMAL
	threadStatus = (hThread != (HANDLE)NULL);
  }
  void stopThread__(){
	  if(threadStatus){
		WaitForSingleObject(hThread, THREAD_TIMEOUT);
		CloseHandle(hThread);
	  }
  }
  
protected:
  IRunnable *_threadObj; 
  static unsigned long __stdcall threadProc(void* ptr) {
    ((IRunnable*)ptr)->run();
    return 0;
  }
  HANDLE hThread;
  int threadStatus;
};


class Model3DAnalysis : IRunnable {
public:
	bool bThreadRunning;
	int nr_files;
	Model3DAnalysis(int nrFiles){
			nr_files = nrFiles;
			pT = new Model3DThread(this);
			if(pT){
				pT->startThread__();
				bThreadRunning = true;
			}
	}

	~Model3DAnalysis(){
		if(pT){
			bThreadRunning = false;
			pT->stopThread__();
		}
	}

	virtual void run();   // this function is redirected to the base class, C++ is marvelous!!:)

protected:
	Model3DThread *pT;

};


#define  FIELD_S 4

struct MIRDATA {
	int x;
	char ho [FIELD_S] ;			// Holzart
	char qk [FIELD_S] ; //übernahme   // Qualität
	char qs [FIELD_S] ;  //sortierung	// Qualität Sortiment
	int la;			// Länge
	int zo;			// Zopf
	int rmi;		// Durchmesser
	int abh;			// Abholzigkeit
	int kr;				// Krümmung
};

class Model3D : IRunnable {
		public:
			MESS* Mess3D;

			MESS Mess3D_Tomolog[1000];   // fix -> this is used for Tomolog 3D, because we have to calc the x,y, not just pointer as for DiShape Data 
			int nSlices; 
			GLVIEW* pGlView;

			CShapeData *pData;	// pointer to the actual object
			CXRayData *pTomo;
			CPhotoData *pPhoto;
			const MACHINECONF *pMC;
			TOMOGEOMETRY GE;

			MIRDATA mirData;
			int mitte;
			int TYPE;
			// Rotating matching fuction stuff
			Model3D *pModelMatch;  // this is the object we're going to match 
			HANDLE hMatchThread;
			DWORD dwIDThread;
			int threadStat;
			bool threadRun;

			// used in the LogInfo3D.c Panel
			LogInfo logInfo_;

			Model3DThread *pThread;
			bool bThreadRunning;

			SOURCE_CALIB sourceCalib; // calibration 3D coords
			KR_LINES lineKr;

			
				 // classes used to store the computed 3D CT slices
			/*
			Image *imageCT;
			 Blob myBlobCT;
			 */
			 unsigned short map[CT_IMG_PIX_SIZE][CT_IMG_PIX_SIZE];	
			 unsigned short mapCT[(CT_IMG_PIX_SIZE/GRID_CELL_SIZE)*(CT_IMG_PIX_SIZE/GRID_CELL_SIZE)];

			 int degMaxKr;
			 



			Model3D(){
				memset(&Mess3D, 0, sizeof(MESS));
				memset(&mirData, 0, sizeof(MIRDATA));
				memset(&lineKr, 0, sizeof(KR_LINES));
				nSlices = 0;
				memset(&pGlView, 0, sizeof(GLVIEW));
				TYPE = 0;
				threadStat = 0;
				bThreadRunning = false;
				LOCK = false;
				degMaxKr = 0;
				lineKr.krAvail = false;
			}
			Model3D(MESS* M3D, int nSli, GLVIEW* pGlV){
				memset(&mirData, 0, sizeof(MIRDATA));
				memset(&lineKr, 0, sizeof(KR_LINES));
				Mess3D = M3D;
				nSlices = nSli;
				pGlView = pGlV;
				TYPE = DISHAPE_WBARK;
				threadStat = 0;
				bThreadRunning = false;
				LOCK = false;
				degMaxKr = 0;
				lineKr.krAvail = false;
			}
			Model3D(CXRayData *pT){
				memset(&mirData, 0, sizeof(MIRDATA));
				memset(&lineKr, 0, sizeof(KR_LINES));
				nSlices = 0;
				pTomo = pT;
				TYPE = XRAY;
				threadStat = 0;
				bThreadRunning = false;
				LOCK = false;
				degMaxKr = 0;
				lineKr.krAvail = false;

				memset(&sourceCalib, 0, sizeof(SOURCE_CALIB));

				// This sets the initial Source coords as calculated by the calibration procedure from Manuel Rezzadore
				sourceCalib.source1.x = -687.1;//-847.1;//-687.1;
				sourceCalib.source1.y = -510.1;//-570.6; //-510.1;
				sourceCalib.source2.x =  707.4;//854.9; //707.4;
				sourceCalib.source2.y = -547.1; //-622.1; //-547.1;

				pMC = pTomo->ExtractConfInfo();
				if(pMC){
					memset(&GE, 0, sizeof(GE));
					GetTomoGeometry(&GE, pMC);
				}
			}

			void initThread(){
				pThread = new Model3DThread(this); 
				if(pThread){
					bThreadRunning = true;
				}
			}

			~Model3D(){
				if(threadStat){
					bThreadRunning = false;
					WaitForSingleObject(hMatchThread, THREAD_TIMEOUT);
					CloseHandle(hMatchThread);
				}
				
			}

			virtual void run(){
					MatchFunc(this);

			}
			int getLength(){
				if(pData != NULL){
					return pData->GetLogLength();
				}else{
					return -1;
				}
			}

			int calcXrayCoords(Model3D *pS, XRAY_DIAMTYPE typeDiam){
				char msg_[256];
				const XRAY_SLICEINFO *info_;
				SLICEINFOCALIB sliceCalib;
				float center_x, center_y, RSx_1, RDx_1, RSx_2, RDx_2, diam1, diam2, stdDiam1, stdDiam2 = 0.0f;
				int  cntAvg = 0;
				memset(&Mess3D_Tomolog, 0, sizeof(Mess3D_Tomolog));

				Mess3D = pS->Mess3D;  // in any case we also set the pointer to DiShape 2 (WOBark)
			
				//fix   we take the nSlices of *pS, 
				nSlices = pS->nSlices;
				center_x = center_y = diam1 = diam2 = 0;
				for(int i = 10; i < nSlices-1; i++){
					if(LOAD_VISUAL_XRAY_CALIB){
						info_ = pTomo->GetSliceInfoByNumber(i);
					}else{
					    info_ = pTomo->GetSliceInfoByPositionRelative(10*(i));  // we extract the relative slice Info
					}
					if(info_ != NULL){
						if(info_->bValid){  // we extract slice info
								memset(&sliceCalib, 0, sizeof(sliceCalib));
								sliceCalib.Sx[0] = (unsigned short)info_->Sx[0];
								sliceCalib.Dx[0] = (unsigned short)info_->Dx[0];
								sliceCalib.Sx[1] = (unsigned short)info_->Sx[1];
								sliceCalib.Dx[1] = (unsigned short)info_->Dx[1];
							 
							FindLogDiametersAndRadialDensityByRezzadore(&sliceCalib, pMC, &sourceCalib);

							center_x = sliceCalib.Xc; //info_->Xc;   //pS->Mess3D[i].zy; //
							center_y =  sliceCalib.Yc;	// pS->Mess3D[i].zy; //
							Mess3D_Tomolog[i].zy = center_y; //pS->Mess3D[i].zy;
							Mess3D_Tomolog[i].zx = center_x; //pS->Mess3D[i].zx;
					//		Mess3D_Tomolog[i].x = pS->Mess3D[i].x;

							if(typeDiam == Optic){
								RSx_1 = info_->RSx[1]*10;  // radii of 1st sensor
								RDx_1 = info_->RDx[1]*10;	
								RSx_2 = info_->RSx[0]*10;  // radii of 2nd sensor
								RDx_2 = info_->RDx[0]*10; 
							}else if(typeDiam == MuMax){
								RSx_1 = info_->RSxMuMax[0]*10;  // radii of 1st sensor
								RDx_1 = info_->RDxMuMax[0]*10;	
								RSx_2 = info_->RSxMuMax[1]*10;  // radii of 2nd sensor
								RDx_2 = info_->RDxMuMax[1]*10; 
							}
							
							diam1 = sliceCalib.D1; //info_->D1;
							diam2 = sliceCalib.D2; //info_->D2;

							stdDiam1 += diam1;
							stdDiam2 += diam2;
							cntAvg++;
			
						//		pS->Mess3D[i].zx = pS->Mess3D[i].zx+center_x;
						//		pS->Mess3D[i].zy = pS->Mess3D[i].zy+center_y;
	
							// fix invert
							double phi = DEGREE_MEASURE_D1;  // was +180
							phi = Grad2radMath(phi);    // fix inverted radii
							Mess3D_Tomolog[i].pointx[DEGREE_MEASURE_D1] = (int)(RSx_1*cos(phi));//+Mess3D_Tomolog[i].zx;  //+(Mess3D_Tomolog[i].zx-center_x)
							Mess3D_Tomolog[i].pointy[DEGREE_MEASURE_D1] = (int)(RSx_1*sin(phi));//+Mess3D_Tomolog[i].zy; //+(Mess3D_Tomolog[i].zy-center_y

							phi = DEGREE_MEASURE_D1+180;
							phi = Grad2radMath(phi);
							Mess3D_Tomolog[i].pointx[DEGREE_MEASURE_D1+180] = (int)(RDx_1*cos(phi));//+Mess3D_Tomolog[i].zx;
							Mess3D_Tomolog[i].pointy[DEGREE_MEASURE_D1+180] = (int)(RDx_1*sin(phi));//+Mess3D_Tomolog[i].zy;

							phi = DEGREE_MEASURE_D2+180;
							phi = Grad2radMath(phi);
							Mess3D_Tomolog[i].pointx[DEGREE_MEASURE_D2+180] = (int)(RSx_2*cos(phi));//+Mess3D_Tomolog[i].zx;
							Mess3D_Tomolog[i].pointy[DEGREE_MEASURE_D2+180] = (int)(RSx_2*sin(phi));//+Mess3D_Tomolog[i].zy;

							phi = DEGREE_MEASURE_D2;    // was + 180
							phi = Grad2radMath(phi);
							Mess3D_Tomolog[i].pointx[DEGREE_MEASURE_D2] = (int)(RDx_2*cos(phi))+Mess3D_Tomolog[i].zx;
							Mess3D_Tomolog[i].pointy[DEGREE_MEASURE_D2] = (int)(RDx_2*sin(phi))+Mess3D_Tomolog[i].zy;


						}else{



						}
					}else{



					}

				}

				stdDiam1 = stdDiam1/cntAvg;
				stdDiam2 = stdDiam2/cntAvg;
				sprintf(msg_, "D1: %f     D2: %f", stdDiam1, stdDiam2); 
				MSG(msg_);
				return 1;
				

			}


			int findBestMatch3D(Model3D *pModelMatch__){
					pModelMatch = pModelMatch__;
					hMatchThread = CreateThread(NULL,
												0,
												(LPTHREAD_START_ROUTINE)MatchFunc,
												this,
												0,
												(LPDWORD)&this->dwIDThread);



					threadStat = (hMatchThread != (HANDLE)NULL);
					SetThreadPriority (hMatchThread, THREAD_PRIORITY_NORMAL);
					return threadStat;
			}

			double CalcDist( MESS* pMess1, MESS* pMess2 ){
				double res = 0.0;
				double resInt = 0;

				if( pMess1 != 0 && pMess2 != 0 ){     
					int* pX1 = &(pMess1->pointx[0]);
					int* pX2 = &(pMess2->pointx[0]);
					int* pY1 = &(pMess1->pointy[0]);
					int* pY2 = &(pMess2->pointy[0]);
					int* pXF = &(pMess1->pointx[360]);
					while( pX1 < pXF )
					{
						resInt += SQR(*pX1 - *pX2) + SQR(*pY1 - *pY2);
						pX1++; pX2++; pY1++; pY2++;
					}
					res = double(resInt)/360.0;
				}
				return res;
			}


			void setLock(bool lock_){
				LOCK = lock_;
			}

			bool getLock(){
				return LOCK;
			}


			void RotateLog(int deg){
				for(int i = 0; i < nSlices; i++){
					RotateSlice(&Mess3D[i], deg);
				}

			}
			

			void RotateSlice( MESS* pMess, int AngleDeg ){
				//this function rotates the slices around the origin.
				//the angle is taken with the goniometric convention.
				if( AngleDeg  == 0 || AngleDeg == 360 )
					return;

				ASSERT( pMess != NULL );   
				ASSERT( AngleDeg >= 0 && AngleDeg < 360 );  //now we support only easy operations

				MESS buf;
				memcpy( &buf, pMess, sizeof(MESS) );
    
				int* pX = &(buf.pointx[0]);
				int* pXF = pX + 360;
				int* pY  = &(buf.pointy[0]);
				double COS = cos( NSP_DegToRad(double(AngleDeg)) );
				double SIN =  sin( NSP_DegToRad(double(AngleDeg)) );
				double tmpX, tmpY;
				int i = AngleDeg;

				while( pX < pXF )
				{
					tmpX = (double)*pX;
					tmpY = (double)*pY;
					pMess->pointx[i] = (int)( COS*tmpX - SIN*tmpY );
					pMess->pointy[i] = (int)( SIN*tmpX + COS*tmpY );
					pX++;
					pY++;
					i++; 
					if( i==360 ) 
						i=0;
				}       
    
				//and now the centre.
				tmpX = (double) pMess->zx;
				tmpY = (double) pMess->zy;
				pMess->zx = (int)( COS*tmpX - SIN*tmpY );
				pMess->zy = (int)( SIN*tmpX + COS*tmpY );
			}



			protected:
				bool LOCK;
		
	};

	/*
		If more complicated constructs are necessary, this class can 
		be extended further  using operator overloading,..some kind
		of abstraction layer;)
  */
	class Vector3D {
	public:
		float x, y, z;
		Vector3D(){

		}
		Vector3D(float x_, float y_, float z_): x(x_), y(y_), z(z_){

		}

		float GetDot(Vector3D &vec){
				return x * vec.x + y * vec.y + z * vec.z;
		}

		Vector3D* GetCross(Vector3D &vec){
				return new Vector3D(y * vec.z - z * vec.y, z * vec.x - x * vec.z, x * vec.y - y * vec.x);
		}

	};


int initRendering(int argc, char **argv);
void OnShowLoadDlg();
void OnLoadMifFiles();
int createTomologObject3D(CXRayData *pTomo);
int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, CShapeData *pShape);
int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, CShapeData *pShape, int index);
int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, int index);
int findMatch(int index);
void setCurrentObj(int index);
void renderLog(int Log_ID); // this function reinitializes the glCallList(id), this 
							// method seems to be faster
int getStatus3D();

  static double Grad2Uhr( double w ){
		// math. Notation -> Uhr
		w = 90 - w;
		while ( w < 0 ) w += 360;
		while ( w > 360 ) w -= 360;
		return w;
	}

	static double Uhr2Grad( double w ){
		// math. Notation -> Uhr
		w = 90 - w;
		while ( w < 0 ) w += 360;
		while ( w > 360 ) w -= 360;
		return w;
	}

	static double radMath2GradUhr( double w ){
		// rad -> Grad
		w = w * 180 / 3.1415926;
		// math. Notation -> Uhr
		return Grad2Uhr( w );
	}

	
	static double radMath2Grad( double w ){
		// rad -> Grad
		w = w * 180 / 3.1415926;
		// math. Notation -> Uhr
		return w;
	}

	static double Grad2radMath(double w){
		w = w * 3.1415926 / 180;
		return w;
	}

	static double gradUhr2RadMath( double w ){
		// Uhr -> math. Notation
		w = Uhr2Grad( w );
		// Grad -> rad
		w = w * 3.1415926 / 180;
		return w;
	}


#endif