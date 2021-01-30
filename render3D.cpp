/**
*	CT (Computed Tomography) Tomolog3DMatching and 3D-Shape Visualization Software (3D-Shape1 measures
*	the log with bark , CT also with bark, and 3D-Shape in debarked status). It provides also an 
*	interface for testing measurment deviations on pre-build testing engines provided by scientific log research.
*	Tested in real-time environment at a sawmill corporation.
*	
*
*
*	
*	Software Author: Armin Costa
*	(c)opyright 2006-2010
*	
*  
*	e-mail: armincosta@hotmail.com
*
*
*	STATEMENT
*	This set of algorithms have been implemented for a custom measurment system in the log industry,
*	Major parts have been also implemented in spare-time activity. In this way, I personally claim some copyright 
*	on this code. 
* 
*	This module implements 3D rendering algorithms, as well as log-data searching-features in order to
*	match heterogenuous log data of all spicies, leading to expanding research possibilities as well as to
*	automatition production businesses.
*.
*	This software renders 3D-Shape Data and CT of the logs with, and without bark. Note that this
*	might not be an easy task as the logs in bark, can eaven have debarked parts.
*	Once the matching has occurred, the XRay projection is also rendered as volume data. In this way
*	We get a true correlation between the 3 different measures. Huge chunks of image files can be
*	analyzed and results being printed to a file, typically comma separated row-alignment.
*
*	file: render3D.cpp
*
*	This source code is distributed under the GNU GENERAL PUBLIC LICENSE Version 3 in the hope it might be helpful
*	
**/
#include "stdafx.h"
#include "render3D.h"
#include "render3DInfoDlg.h"
#include "..\lib\SLXData\shapedata.h"
#include "..\lib\SLXData\MirData.h"
#include "..\lib\SLXData\XRayStructs.h"
#include "..\lib\3dstructs.h"
#include "..\Interface\loadObjdlg.h"
#include "..\Common\Interface.h"
#include "..\Common\parametri.h"
#include "..\Core\TomologUtil.h"

//#include "..\lib\DiShape3d.h"  // this is just the newer version of DiShape3D structs, for reading DiShape exported data
#include "..\Core\osdir.h"

#include "..\lib\render3DTomolog.h"   // this is the magic interface
/*
#include <Magick++.h>  	// we use Image Libs from ImageMagick to operate on image formats
#include <stdio.h>
			
using namespace Magick;

*/

using namespace std;



/*  3D Model params  */
GLfloat fovy=15;
GLfloat near_plane=1; // 2
GLfloat far_plane=1000; // 10
GLfloat alfa=28;
GLfloat beta=10;
GLfloat z_angle = 0.1f;
GLfloat radius;
GLfloat Zoom;

static double theta = 61.4;

#define DRAFT  0
#define MEDIUM 1
#define BEST   2
int drawquality = MEDIUM;
int spincamera = TRUE;
int cameradirection = 0;
double updownrotate = -78;
//#########################################################################################

DATAOBJECT LoadMifObj[MAX_LOGS3D];   // this is used to Load a 3D model from a mif file
// Log Info
int getLogDiam_WBark(int step,  int deg_measure);
int getLogDiam_WOBark(int step,  int deg_measure);
int getLogDiam_WBarkTOMOLOG(int step, int deg_measure);
int getLogDiam_WOBarkTOMOLOG(int step, int deg_measure, XRAY_DIAMTYPE alg);

const XRAY_SLICEINFO* getXRAY_SLICEINFO(int step);
int getStdLogDiam_WBark(int step, int range, int deg_measure,  bool cut_off);
int getStdLogDiam_WOBark(int step, int range, int deg_measure,  bool cut_off);
int getStdLogDiam_WBarkTOMOLOG(int step, int range, int deg_measure,  bool cut_off);
int getStdLogDiam_WOBarkTOMOLOG(int step, int range, int deg_measure,  XRAY_DIAMTYPE alg, bool cut_off);

int getLengthWBark();
int getLengthWOBark();

int getMinimumDiamTOMO(int start, int end, int &off_ret, XRAY_DIAMTYPE type);  // returns the min avg Diam (M1+M2) and sets the len offset found
int getMinimumDiamDiShapeWBark(int start, int end, int &off_ret);
int getMinimumDiamDiShapeWOBark(int start, int end, int &off_ret);

void showLogData(bool writeToFile);  // this function writes results to file
 
void populateONORM(HFAData &data);

//###########################################
// This macros are used in the function that runs the Peintinger Formula
#define FAKTOR_FICHTE		95931
#define OFFSET_FICHTE		-256182
#define FAKTOR_LAERCHE		88260
#define OFFSET_LAERCHE		-208500
#define FAKTOR_KIEFER       95075L
#define OFFSET_KIEFER       -173823L

int getPeintingerMeasure(int diam, char type);
//###########################################

void initGL();

GLvoid drawTomologLog();
GLvoid drawWireframeLog();
GLvoid drawWireframeLogMatch();
GLvoid LogDiskKr(int Slice_, int k, float alpha,  MESS *s1,MESS *s2, float mitte, GLVIEW* pGlView, DWORD Type_);
GLvoid LogDisk1(int Slice_, int k, float alpha,  MESS *s1,MESS *s2, float mitte, GLVIEW* pGlView, DWORD Type_);
GLvoid LogDisk2(int Slice_, int k, float alpha,  MESS *s1,MESS *s2, float mitte, GLVIEW* pGlView, DWORD Type_);
GLvoid LogDiskTomo(int Slice_, int k, float cnt, MESS *s1,MESS *s2, float mitte,GLVIEW* pGlView, DWORD Type_);
GLvoid Disk(int k,MESS *s);
GLvoid GetNormal(GLfloat v1[3],GLfloat v2[3],GLfloat n[3]);
void normCrossProd(float v1[3], float v2[3], float out[3]);
void normalize(float v[3]);
void Display(void);
void CreateEnvironment(void);
void MakeGeometry(void);
void MakeLighting(void);
void MakeCamera(int,int,int);
void MakeMyFont();
void renderLog(int Log_ID);
void HandleKeyboard(unsigned char key,int x, int y);
void HandleSpecialKeyboard(int key,int x, int y);
void HandleMouse(int,int,int,int);
void HandleMainMenu(int);
void HandleSpeedMenu(int);
void HandleVisibility(int vis);
void HandleIdle(void);
void DrawTextXY(double,double,double,double,char *);
void GiveUsage(char *);
void printString(char *s);
void myPrint(float x, float y, float z, char* format, ...);
GLvoid draw3DAxis();

// 3D Tomolog stuff
GLvoid exportTomologVolume();

BOOL ReadMessung(char *path, LOG *plog);  // this function servers to read 3D DiShape Data Files as exported in DiShape

int compareMin(const void * a, const void * b);  // this is the Callback function for qsort, used to extract minimums
// defined in shapedata.h
static double CalcDist( MESS* pMess1, MESS* pMess2 );
int UP = 0;
#define TRUE  1
#define FALSE 0
#define PI 3.141592653589793238462643
/*
	this Data types are used to load *.mif files
  */
	CShapeData* p3D1 = NULL;
   CShapeData* p3D2 = NULL;
   CXRayData* pXRay = NULL;
    CXRayData *pXR;   // we use this to load raw tomolog obj (.mif)
   CStringList *m_pFileToLoadList;
   CString mif_file;
   const char *pMif;
   SYSTEMTIME myTime;
   POSITION m_CurrentPosition;
    //load the given file
    CLoadObjDlg LoadObjDlg;

	HFAData dataSet;   // we use this struct for evaluation purposes
	KRData  dataSetKr;

/*
	these are the main objects 
  */ 
////////////////////////////////
	 Model3D *myLogs[MAX_LOGS3D*2]; // because here we store 3D1 and 3D2 (DiShape)
	 int logCnt = 0;

	 Model3D *myXRayLog[MAX_LOGS3D]; // no data in here, just pointers
	 int logCntXRay = 0;

	 int logCntXray = 0;
	 int currentLog = 0;
	 int logMatch = 0;

	 TomoRenderer *ctRender;

///////////////////////////////////
	 bool defaultVModel = true;
	 bool frontVModelP = false;
	
	 bool showWBark = true;	// 7
	 bool showWOBark = true;	//8
	 bool showXRay = false;

	 bool showOnlyMark = false;
	 bool DlgLogData = false;
	 int INFO_WIN_ID = 0;
	 CShapeData tmpShape; // I work on a copy to keep constistet
	 MESS tmpMess;

	 double areaStart = 0;
	 double areaEnd = areaStart+SAMPLE_INTERVAL;

	 CLogInfo3D logInfoDlg;

	 Model3DAnalysis *pModelAnalysis;
	 static CString fileList[MAX_MIF_FILES];
	 static CString showfileList[MAX_MIF_FILES];

GLvoid *font_style = GLUT_BITMAP_TIMES_ROMAN_24;

static bool LOCKED = false;

void myPrint(float x, float y, float z, char* format, ...);


void exportLogDataKr(bool writeToFile);  // this is just another Log Analysis for  curvature

BOOL getKruemmungMittellinieBx(LOG *log, int off_start, int off_end);
BOOL getKruemmungOberflaecheLinz(LOG *log, int off_start, int off_end, int degMaxKr, int diam__); // Microtec Linz
BOOL getKruemmungMittellinieLinz(LOG *log, int off_start, int off_end, int diam__, double *dWinkel);
BOOL getKrAngle(LOG *log, int off_start, int off_end, int diam__, double *dWinkel);
void smoothLog(Model3D *pModel);

void VisualizeKruemmung();


bool isMeasureOverrun(const XRAY_SLICEINFO *pS);


int hasEinseitigerWurzelAn(LOG *log);
int	BestimmeIndex(LOG *log,int lv,int lb,int *iv,int *ib);
LOG myL;   // Strukture used to calc the Curvature  -> just because we want to recycle code:)


int createTomologObject3D(CXRayData *pTomo){
	if(logCntXRay < MAX_LOGS3D){
		myXRayLog[logCntXray] = new Model3D(pTomo); 
		myXRayLog[logCntXray]->mitte = 0;
		 logCntXray++;
		return 1;
	}else{
		return 0;	
	}
}

int createTomologObject3D(CXRayData *pTomo,  int nSlices, int index){
	if(index < MAX_LOGS3D){
		myXRayLog[index] = new Model3D(pTomo); 
		myXRayLog[index]->nSlices = nSlices;
		myXRayLog[index]->mitte = nSlices/2;

//		myXRayLog[index]->pTomo = pTomo;
		if(myXRayLog[0] != NULL){
			myXRayLog[index]->calcXrayCoords(myXRayLog[0], (XRAY_DIAMTYPE)0);   //we take DiShape 1 as reference coords
		}
		logCntXray++;
		return 1;
	}else{
		return 0;	
	}
}

Model3D *getXRayLog(){
		return (myXRayLog[currentLog] != NULL) ? myXRayLog[currentLog] : myXRayLog[currentLog] = new Model3D();

}

int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, CShapeData *pShape){
	if(logCnt < MAX_LOGS3D*2){
		myLogs[logCnt] = new Model3D(Mess3D, nSlices, pGlView); 
		myLogs[logCnt]->mitte = nSlices/2;
		myLogs[logCnt]->pData = pShape; logCnt++;
		return 1;
	}else{
		return 0;
	}	
}

int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, CShapeData *pShape, int index){
	if(logCnt < MAX_LOGS3D*2){
		myLogs[index] = new Model3D(Mess3D, nSlices, pGlView); 
		myLogs[index]->mitte = nSlices/2;
		myLogs[index]->pData = pShape; 
		return 1;
	}else{
		return 0;
	}	
}

int createObj3D( MESS* Mess3D, int nSlices, GLVIEW* pGlView, int index ){
	if(index < MAX_LOGS3D*2){
		myLogs[index] = new Model3D(Mess3D, nSlices, pGlView); 
		return 1;
	}else{
		return 0;
	}
}

void setCurrentObj(int index){
	if(index < MAX_LOGS3D*2){
		currentLog = index;
	}
}

int getStatus3D(){
	return UP;
}

/*
	This function takes care of reinitializing the proper rendering glCallList 
	Param: 1 == Log W Bark
		   2 == Log WO Bark
  */
void renderLog(int Log_ID){
	if(!COLLECT_DATA){ // to ensure that this is not called in Loading Operation (we want to load && analyze 10000 objs) 
		if(Log_ID == 1){
			glNewList(Log_ID, GL_COMPILE);
			drawWireframeLog();
			glEndList();
		}else if(Log_ID == 2){
			glNewList(Log_ID, GL_COMPILE);
			drawWireframeLogMatch();
			glEndList();
		}else if(Log_ID == 3){
			glNewList(Log_ID, GL_COMPILE);
			drawTomologLog();
			glEndList();

		}else if(Log_ID == 100){
			glNewList(1, GL_COMPILE);
			drawWireframeLog();
			glEndList();
			glNewList(2, GL_COMPILE);
			drawWireframeLogMatch();
			glEndList();
			glNewList(3, GL_COMPILE);
			drawTomologLog();
			glEndList();
			
		}
	}

}


int findMatch(int index){
	int i, x,y, deg, len, pInd, p, i_p;
	int s_dist = 0;
	int min_dist = 100000000;
	int min_dist2 = 100000000;
	double min_[MAX_LOGS3D];
	char msg[100];
	double avgDist;

	if(myLogs[index] != NULL){
		len = myLogs[index]->nSlices/2;

		if(TRUST_SLX_MATCHING){
			i = index+1; // we take the Log matched by SLX
			if(myLogs[i] != NULL && (myLogs[i]->TYPE != XRAY)){
				if(i != index){ // we exclude the exact object	
					myLogs[i]->pData->AlignToObject( myLogs[index]->pData, len); // we align the obj
					s_dist = (int)CalcDist(&myLogs[i]->Mess3D[len], &myLogs[index]->Mess3D[len]);
					min_dist = min(min_dist, s_dist);
					if(min_dist == s_dist){
						pInd = i;
						logMatch = pInd;
					}
				}
			}
		}else{
			for(i = 0; i < MAX_LOGS3D*2; i++){
				if(myLogs[i] != NULL && (myLogs[i]->TYPE != XRAY)){
					if(i != index){ // we exclude the exact object	
						myLogs[i]->pData->AlignToObject( myLogs[index]->pData, len); // we align the obj
						s_dist = (int)CalcDist(&myLogs[i]->Mess3D[len], &myLogs[index]->Mess3D[len]);
						min_dist = min(min_dist, s_dist);
						if(min_dist == s_dist){
							pInd = i;
							logMatch = pInd;
						}
					}
				}
			}
		}
		if(!COLLECT_DATA){
			sprintf(msg, "min Dist %d | index:%d ", min_dist, pInd);
			MSG(msg);
		}
	}
	return pInd;
}

void getDataDiameters(){
/*
        int len, Min0, Min90, Max0, Max90;

        len = myLogs[currentLog]->getLength();

        //first end:
        GetMinMaxDiameterInRange( pDmAnalysisData->nTopSearchStart, 
            pDmAnalysisData->nTopSearchStart + pDmAnalysisData->nTopLengRange,
            0,
            &Min0, &Max0 );
        GetMinMaxDiameterInRange( pDmAnalysisData->nTopSearchStart, 
            pDmAnalysisData->nTopSearchStart + pDmAnalysisData->nTopLengRange,
            90,
            &Min90, &Max90 );
        MinDiam1 = min( Min0, Min90 ); MaxDiam1 = max( Max0, Max90 ) ;	
		*/

}


int getLogDiam_WBark(int step, int deg_measure){
	int d1;
	if(myLogs[currentLog] != NULL){
		d1 = (myLogs[currentLog]->pData)->GetDiameterAt(step, deg_measure, bClambDiShapeDiam);//((myLogs[currentLog]->Mess3D[4]), 40, false);
		return d1;
	}else{
		return 0;
	}
}

int getLogDiam_WOBark(int step, int deg_measure){
	int d2;
	if(myLogs[logMatch] != NULL){
		d2 =  (myLogs[logMatch]->pData)->GetDiameterAt(step, deg_measure, bClambDiShapeDiam);
		return d2;;
	}else{
		return 0;
	}	
}

int getLogDiam_WBarkTOMOLOG(int step, int deg_measure){
	if(myXRayLog[currentLog] != NULL){
		if(deg_measure == DEGREE_MEASURE_D1){
			if(USE_OFFSET){
				return (int)((myXRayLog[currentLog]->pTomo)->GetDiameter (step, 1, (XRAY_DIAMTYPE)0)*10)+TOMOLOG_WB;
			}else{
				return (int)(myXRayLog[currentLog]->pTomo)->GetDiameter (step, 1, (XRAY_DIAMTYPE)0)*10;
			}
		}else if(deg_measure == DEGREE_MEASURE_D2){
			if(USE_OFFSET){
			   return (int)((myXRayLog[currentLog]->pTomo)->GetDiameter (step, 0, (XRAY_DIAMTYPE)0)*10)+TOMOLOG_WB2;
			}else{
			   return (int)(myXRayLog[currentLog]->pTomo)->GetDiameter (step, 0, (XRAY_DIAMTYPE)0)*10;
			}
		}
	}else{
		return 0;	
	}
}

int getLogDiam_WOBarkTOMOLOG(int step, int deg_measure, XRAY_DIAMTYPE alg){
	if(myXRayLog[currentLog] != NULL){
		if(deg_measure == DEGREE_MEASURE_D1){
			if(USE_OFFSET){
				return (int)((myXRayLog[currentLog]->pTomo)->GetDiameter (step, 0, alg)*10)+TOMOLOG_WOB;
			}else{
				return (int)(myXRayLog[currentLog]->pTomo)->GetDiameter (step, 0, alg)*10;
			}
		}else if(deg_measure == DEGREE_MEASURE_D2){
			if(USE_OFFSET){
				return (int)((myXRayLog[currentLog]->pTomo)->GetDiameter (step, 1, alg)*10)+TOMOLOG_WOB2;
			}else{
				return (int)(myXRayLog[currentLog]->pTomo)->GetDiameter (step, 1, alg)*10;
			}
		}
	}else{
		return 0;	
	}
}



int getStdLogDiam_WBark(int step, int range, int deg_measure, bool cut_off){
	int  median, acc_m, start_off, end_off, i, cnt, off;
	start_off = step-range;
	end_off = step+range;
	acc_m = 0;
	cnt = 0;
	int diam;
	int size = range*2+1; 
	static int val[1000];
	int s_off, e_off;

	if(!cut_off){
		if(myLogs[currentLog] != NULL){
			for(i = start_off; i <= end_off; i++){
				acc_m += getLogDiam_WBark(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;
		}else{
			return 0;
		}
	}else{
		if(myLogs[currentLog] != NULL){
			off = 0;
			for(i = start_off; i <= end_off; i++){ 
								diam =  (int)getLogDiam_WBark(i, deg_measure);
								val[off] = diam;     off++;
			}
			qsort(val, size-1, sizeof(int), compareMin);
			// we exclude the 25% of min. and 25% of max.
			s_off = (int)(size-1/2)/2;
			e_off = size-s_off;
			for(i = s_off; i <= e_off; i++){
				acc_m += getLogDiam_WBark(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;
		}
	}
}

int getStdLogDiam_WOBark(int step, int range, int deg_measure, bool cut_off){
	int  median, acc_m, start_off, end_off, i, cnt, off;
	start_off = step-range;
	end_off = step+range;
	acc_m = 0;
	cnt = 0;
	int diam;
	int size = range*2+1; 
	static int val[1000];
	int s_off, e_off;

	if(!cut_off){
		if(myLogs[logMatch] != NULL){
			for(i = start_off; i <= end_off; i++){
				acc_m += getLogDiam_WOBark(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;
		}else{
			return 0;
		}
	}else{
		if(myLogs[logMatch] != NULL){
			off = 0;
			for(i = start_off; i <= end_off; i++){ 
								diam =  (int)getLogDiam_WOBark(i, deg_measure);
								val[off] = diam;     off++;
			}
			qsort(val, size-1, sizeof(int), compareMin);
			// we exclude the 25% of min. and 25% of max.
			s_off = (int)(size-1/2)/2;
			e_off = size-s_off;
			for(i = s_off; i <= e_off; i++){
				acc_m += getLogDiam_WOBark(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;

		}
	}

}

int getStdLogDiam_WBarkTOMOLOG(int step, int range, int deg_measure,  bool cut_off){
	int  median, acc_m, start_off, end_off, i, cnt, off;
	start_off = step-range;
	end_off = step+range;
	acc_m = 0;
	cnt = 0;
	int diam;
	int size = range*2+1; 
	static int val[1000];
	int s_off, e_off;

	if(!cut_off){
			for(i = start_off; i <= end_off; i++){
				acc_m += getLogDiam_WBarkTOMOLOG(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;
	}else{
			off = 0;
			for(i = start_off; i <= end_off; i++){ 
								diam =  (int)getLogDiam_WBarkTOMOLOG(i, deg_measure);
								val[off] = diam;     off++;
			}
			qsort(val, size-1, sizeof(int), compareMin);	
			// we exclude the 25% of min. and 25% of max.
			s_off = (int)(size-1/2)/2;
			e_off = size-s_off;
			for(i = s_off; i <= e_off; i++){
				acc_m += getLogDiam_WBarkTOMOLOG(i, deg_measure);  cnt++;
			}
			median = acc_m/cnt;
			return median;
	}
}


int getStdLogDiam_WOBarkTOMOLOG(int step, int range, int deg_measure,  XRAY_DIAMTYPE alg, bool cut_off){
	int  median, acc_m, start_off, end_off, i, cnt, off;
	start_off = step-range;
	end_off = step+range;
	acc_m = 0;
	cnt = 0;
	int diam;
	int size = range*2+1; 
	static int val[1000];
	int s_off, e_off;

	if(!cut_off){
			for(i = start_off; i <= end_off; i++){
				acc_m += getLogDiam_WOBarkTOMOLOG(i, deg_measure, alg);  cnt++;
			}
			median = acc_m/cnt;
			return median;

	}else{
			off = 0;
			for(i = start_off; i <= end_off; i++){ 
								diam =  (int)getLogDiam_WOBarkTOMOLOG(i, deg_measure, alg);
								val[off] = diam;     off++;
			}
			qsort(val, size-1, sizeof(int), compareMin);

			
			// we exclude the 25% of min. and 25% of max.
			s_off = (int)(size-1/2)/2;
			e_off = size-s_off;

			for(i = s_off; i <= e_off; i++){
				acc_m += getLogDiam_WOBarkTOMOLOG(i, deg_measure, alg);  cnt++;
			}
			median = acc_m/cnt;
			return median;
	}

}


const XRAY_SLICEINFO* getXRAY_SLICEINFO(int step){
	static const XRAY_SLICEINFO *pInfo; 
	Model3D *pM = getXRayLog();
	pInfo = (pM->pTomo)->GetSliceInfoByPositionRelative(step);  // because we take mm
	return pInfo;
}

int getLengthWBark(){
	if(myLogs[currentLog] != NULL){
		return myLogs[currentLog]->getLength();
	}else{
		return 0;
	}

}

int getLengthWOBark(){
	if(myLogs[logMatch] != NULL){
		return myLogs[logMatch]->getLength();
	}else{
		return 0;
	}
}

GLvoid drawLogInfo(int step){
		glColor3f(0.0, 0.0, 3.0);
	int r1 = getLogDiam_WBark(step, DEGREE_MEASURE_D1);
	int r2 = getLogDiam_WOBark(step, DEGREE_MEASURE_D1);
	myPrint(1.0, 1.0, 0, "Diam with Bark: %d\n", r1);
	myPrint(1.0, 1.0, 0, "Diam without Bark: %d\n", r2);	
}

 
static double CalcDist( MESS* pMess1, MESS* pMess2 ){
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


int initRendering(int argc, char **argv)
{
   int i,j,depth;
   int mainmenu,speedmenu;
   unsigned int subWin;
   drawquality = MEDIUM;

   glutInit(&argc,argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

   glutCreateWindow("MiCROTEC.at  SLX-Tomolog    Author: Armin Costa");
   glutReshapeWindow(550, 500);
   glutPositionWindow(100, 100);
   glutDisplayFunc(Display);
   glutVisibilityFunc(HandleVisibility);
   glutKeyboardFunc(HandleKeyboard);
   glutSpecialFunc(HandleSpecialKeyboard);
   glutMouseFunc(HandleMouse);

  // glutCreateWindow(GetActiveWindow(), "MiCROTEC.at  -- Tomolog");
  // glutCreateSubWindow(NULL, 100, 100, 500, 500);

   CreateEnvironment();
   speedmenu = glutCreateMenu(HandleSpeedMenu);
   glutAddMenuEntry("Slow",1);
   glutAddMenuEntry("Medium",2);
   glutAddMenuEntry("fast",3);

   mainmenu = glutCreateMenu(HandleMainMenu);
   glutAddMenuEntry("Toggle camera spin",1);
   glutAddMenuEntry("Toggle ball bounce",2);
   glutAddSubMenu("Ball speed",speedmenu);

   glutAddMenuEntry("3D Mode 1",4);
   glutAddMenuEntry("3D Mode 2",5);
   glutAddMenuEntry("fit Model", 6);
   glutAddMenuEntry("fit Exact Model", 20);
   glutAddMenuEntry("w Bark", 7);
   glutAddMenuEntry("wo Bark", 8);
   glutAddMenuEntry("export X-Ray 3D data", 13);
   glutAddMenuEntry("show both", 9);
   glutAddMenuEntry("show only Mark WO Bark", 40);
   glutAddMenuEntry("show CT", 41);
   glutAddMenuEntry("show Data", 10);
   glutAddMenuEntry("Geometry", 50);
   glutAddMenuEntry("Load..", 30);
   glutAddMenuEntry("Krümmung...", 66);
   glutAddMenuEntry("Quit",100);
   glutAttachMenu(GLUT_RIGHT_BUTTON);

	initGL();

	// we swap the order in order to avoid overlays
	renderLog(1);  // with bark
	renderLog(2);	// without bark, used to match
	renderLog(3); // this is the visual approximation of the Tomolog contour

		MSG("UP");
   glutMainLoop();
   UP = 1;



   ctRender = new TomoRenderer();  // this is the 3D CT rendering module
   
   return(0);
}



void initGL(){
  GLfloat lightPosition[] = {-3.0, 0.0, 1.5, 1.0};
    GLfloat matSpecular[] = {0.8, 0.8, 0.8, 1.0};
    GLfloat matShininess[] = {80.0};

	 int width = 500;
	 int height = 500;
 	 width-=20;
 	 height-=20;
	 Zoom = 20.0f;
    glViewport( 10, 10, width, height ); 

    glClearColor(0.9f, 0.9f, 1.0f, 1.0f) ;
	glShadeModel(GL_SMOOTH);	
	glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
    glMaterialfv(GL_FRONT, GL_SHININESS, matShininess);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_LIGHTING) ;
    glEnable(GL_LIGHT0) ;

    glMatrixMode( GL_PROJECTION );     
    gluPerspective( 45.0, (GLfloat)(width / height), 1.0, 12.0 ); 
    glMatrixMode( GL_MODELVIEW ); 


}


/*
	this function exports a set of projection slices into a series of renderable files.
	The Algorithm constructs from two single 2D density profiles (Sensor 1 && Sensor 2) a 2D Image which incorporates
	a true 3D model approximation with only 2D trigonimetric measures. This model can be adapted to up to N Sensors,
	but also just 3 or 4 Sensors would reveal into a cheap and fine 3D CT.
	This algorithm can be extended and optimized futher by implementing a sort of cell based avaraged Grid.
	This function also allows to export RAW 16bit image Data.
	NOTE: For efficiency reasons a lower resolution Grid has been constructed, by takig the avg of every Cell
	This Algorithm has been designed and invented by Armin Costa
	
  */
class GridCell {
public:
	GridCell(){
		width = 0;
		height = 0;
		valuePix = 0;
		memset(&rowAvg, 0, sizeof(rowAvg));
		int avg = 0;
		int counter = 0;
		int rowCounter = 0;
	}
	DWORD width;
	DWORD height;
	unsigned short valuePix;
	unsigned short rowAvg[GRID_CELL_SIZE*GRID_CELL_SIZE];
	int avg;
	int counter;
	int rowCounter;

};
const int gridSize = CT_IMG_PIX_SIZE/GRID_CELL_SIZE;
GridCell myGrid[gridSize*gridSize];



BOOL SimpleInPlaceSmooth( float* samples, int nLength, int nOrder )
{
    if( nLength > 0 && nLength > nOrder )
    {
        float *pKernel  = nspsMalloc( nOrder );
        float *pSupport = nspsMalloc( (nLength + nOrder -1 ) );
        
        //prepare the kernel and the output...
        nspsbSet( 1.0f / (float)nOrder, pKernel, nOrder );
        nspsbZero( pSupport , nLength + nOrder -1 );
        
        //convolve the signal and the kernel...
        nspsConv( samples, nLength, pKernel, nOrder, pSupport );
		
        
        //copy the smoothed data into the original input array.
        nspsbCopy( pSupport + (nOrder / 2), samples, nLength );
        
        nspFree( pKernel );
        nspFree( pSupport );
        return TRUE;        
    }
    return FALSE;    
}

GLvoid exportTomologVolume(){
	static const XRAY_SLICEINFO *pInfo; 
	int logLen = 0;
	int sizeF = 0;
	const float *pProj_1;
	const float *pProj_2;
	const size_t size = CT_IMG_PIX_SIZE*CT_IMG_PIX_SIZE;//768*768;
	unsigned short valuePix_1[CT_IMG_PIX_SIZE];
	unsigned short valuePix_2[CT_IMG_PIX_SIZE];
	unsigned short tmp = 0;
	char path_[256];
	char buff[20];
	const XRAY_HEADER *pHead;
	FILE *pF;

	std::string formatS = "PGM";  // this is PGM

	if(myXRayLog[currentLog] != NULL){
		memset(&valuePix_1, 0, sizeof(valuePix_1));
		memset(&valuePix_2, 0, sizeof(valuePix_2));
		pHead = (myXRayLog[currentLog]->pTomo)->ExtractHeaderInfo();
	    DWORD   TotalLines    = pHead->NumOfRows;
		int     ByteLineWidth = pHead->RowByteSize;
		int     NumPixel      = pHead->NumPixels;
		logLen = (myXRayLog[currentLog]->pTomo)->GetLogLength();
		
		pProj_1 = (myXRayLog[currentLog]->pTomo)->GetProjection(0);
		pProj_2 = (myXRayLog[currentLog]->pTomo)->GetProjection(1);


		int sliceCnt = 1;
		int sStart =  TotalLines/2;
		int sEnd =  TotalLines;

		sStart = sStart;
	//	sEnd =  sStart+200;

		myGrid[0].height = GRID_CELL_SIZE;
		myGrid[0].width = GRID_CELL_SIZE;
		
		
		pInfo = (myXRayLog[currentLog]->pTomo)->GetSliceInfoByPositionAbsolute(1500);  // because we take mm
		
		// if we want to sample raw Tomolog signals, we have to use this routine
	//	pInfo = myXRayLog[currentLog]->pTomo->GetSliceInfoByNumber(sEnd/2);
		//	int len = myXRayLog[currentLog]->getLength();
		
		// we have to find the Pixel boudaries of the log -> we have to optimize 
		int bound_s_1, bound_d_1, bound_s_2, bound_d_2 = 0;
		bound_s_1 = pInfo->Sx[0];
		bound_d_1 = pInfo->Dx[0];

		bound_s_2 = pInfo->Sx[1];
		bound_d_2 = pInfo->Dx[1];

	//	  --> this seems to work better because we actually see an irregular countours
		// we add an offset to the bounds, just to capture the whole log
		bound_s_1 = 0; //bound_s_1-20;
		bound_d_1 = 0; //bound_d_1+20;
		bound_s_2 = 0;
		bound_d_2 = 0; //bound_d_2+20;
		

		for(int x = sStart; x < sEnd; x++){
				sprintf(path_, "%s\\proj.%d", CT_FILES, sliceCnt);   sliceCnt++;
				string ctSlice(path_);
				for(int i = 0; i < NumPixel; i++){
					if((pProj_1+(x*NumPixel)+i) != NULL && (pProj_2+(x*NumPixel)+i) != NULL){
						valuePix_1[i] = (unsigned short)((*(pProj_1+(x*NumPixel)+i))*100);
						valuePix_2[i] = (unsigned short)((*(pProj_2+(x*NumPixel)+i))*100);	
						/*
							if(valuePix_1[i] > 65536 || valuePix_1[i] > BROKEN_PIX_THRES){
								valuePix_1[i] = 0; //we set the value low
							}
							if(valuePix_2[i] > 65536 || valuePix_2[i] > BROKEN_PIX_THRES){
								valuePix_2[i] = 0; //we set the value low
							}
							*/
					}
					
				}
				// NOTE: 
				if(x == sEnd/2){
					printf("mega");
				}
				// lets build a 2D image from both projections  -> NOTE: this is just an apporzimation
				memset(&myXRayLog[currentLog]->map, 0, sizeof(myXRayLog[currentLog]->map));
				memset(&myXRayLog[currentLog]->mapCT, 0, sizeof(myXRayLog[currentLog]->mapCT));

				// we apply the algorithm to the bound interval
				for(int z = bound_s_1; z < NumPixel-bound_d_2;z++){
					for(int u = bound_s_2; u < NumPixel-bound_d_2;u++){
							tmp = (unsigned short)((valuePix_1[z]+valuePix_2[u])/2);  // this is the avg
							myXRayLog[currentLog]->map[z][u] = tmp;	
						//	myXRayLog[currentLog]->mapCT[u+z*NumPixel] = myXRayLog[currentLog]->map[z][u];
					    	myXRayLog[currentLog]->mapCT[u*NumPixel+z] = myXRayLog[currentLog]->map[z][u];
							
					}
				}


				//SimpleInPlaceSmooth(
	

				
				pF = fopen(path_, "wb");
				if(pF){
					// we just write RAW gray mapped data
					fwrite((void*)(&myXRayLog[currentLog]->mapCT), 1, sizeof(myXRayLog[currentLog]->mapCT), pF);
					fflush(pF);
					Sleep(200);
					fclose(pF);
					
				}
				
					

		}

	}
}

GLvoid drawTomoGeom(int step){
	static char m_[256];
	const MACHINECONF *pMC = NULL;
	const XRAY_SLICEINFO *sInfo = NULL;
	XRAY_SLICEINFO tmpS;
	SLICEINFOCALIB sliceCalib;
	TOMOGEOMETRY GE;
	float x, y, det_x1, det_y1, det_x2, det_y2, proj_x1, proj_y1, proj_x2, proj_y2, center_Mx, center_My;


	pMC = myXRayLog[currentLog]->pMC;
	memcpy(&GE, &myXRayLog[currentLog]->GE, sizeof(GE)); 

	memset(&tmpS, 0, sizeof(tmpS));
	memset(&sliceCalib, 0, sizeof(sliceCalib));

	if(LOAD_VISUAL_XRAY_CALIB){
			sInfo = myXRayLog[currentLog]->pTomo->GetSliceInfoByNumber(step);
	}else{
		sInfo = myXRayLog[currentLog]->pTomo->GetSliceInfoByPositionRelative(step*10);
	}
		if(sInfo != NULL){
			//	memcpy(&tmpS, sInfo, sizeof(sInfo));
			//	FindSliceCenterAndTangent(&tmpS, &GE, pMC);
			    sliceCalib.Sx[0] = sInfo->Sx[0];
				sliceCalib.Dx[0] = sInfo->Dx[0];
				sliceCalib.Sx[1] = sInfo->Sx[1];
				sliceCalib.Dx[1] = sInfo->Dx[1];
				 
				FindLogDiametersAndRadialDensityByRezzadore(&sliceCalib, pMC, &(myXRayLog[currentLog]->sourceCalib));
			//	MSG("end");


				// Source 1 and Detector 1
				center_Mx = (float)GE.Cx[0].x/1000;
				center_My = (float)GE.Cx[0].y/1000;

				x = (float)sliceCalib.S[1].x/1000; //GE.S[0].x/1000; //pMC->S[0].x;
				y = (float)sliceCalib.S[1].y/1000; //GE.S[0].y/1000; //pMC->S[0].y;
				det_x1 = (float)GE.DetFirstPix[0].x/1000;
				det_y1 = (float)GE.DetFirstPix[0].y/1000;
				det_x2 = (float)GE.DetLastPix[0].x/1000;
				det_y2 =  (float)GE.DetLastPix[0].y/1000;

				proj_x1 = x,
				proj_y1 = y;
				sprintf(m_, "Q1:   %f, %f", x, y);
		//		MSG(m_);
				proj_x2 = (float)sliceCalib.Left[1].x/1000; //sInfo->Left[0].x/1000;
				proj_y2 =  (float)sliceCalib.Left[1].y/1000; //sInfo->Left[0].y/1000;

							glColor3f(0.0, 0.0, 0.5);
							glBegin(GL_POLYGON);   // Source 1
							glVertex3f(x,y,0.0);
							glVertex3f(x+0.05,y,0.0);
							glVertex3f(x+0.05,y-0.05,0.0);
							glVertex3f(x,y-0.05,0.0);
							glEnd();

							glBegin(GL_LINE_STRIP);	// Detector 1
							glVertex3f(det_x1, det_y1, 0.0);
							glVertex3f(det_x2, det_y2, 0.0);
							glEnd();

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							proj_x2 = (float)sliceCalib.Right[1].x/1000; //sInfo->Right[0].x/1000;
							proj_y2 =  (float)sliceCalib.Right[1].y/1000;//sInfo->Right[0].y/1000;

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							proj_x2 = (float)sliceCalib.Ctr[1].x/1000; //sInfo->Right[0].x/1000;
							proj_y2 =  (float)sliceCalib.Ctr[1].y/1000;//sInfo->Right[0].y/1000;

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							glBegin(GL_POLYGON);   //Center
							glVertex3f(center_Mx,center_My,0.0);
							glVertex3f(center_Mx+0.05,center_My,0.0);
							glVertex3f(center_Mx+0.05,center_My-0.05,0.0);
							glVertex3f(center_Mx,center_My-0.05,0.0);
							glEnd();

				
				// Source 2 and Detector 2
				center_Mx = (float)GE.Cx[1].x/1000;
				center_My = (float)GE.Cx[1].y/1000;
				x = (float)sliceCalib.S[0].x/1000; //GE.S[1].x/1000; //pMC->S[1].x;
				y = (float)sliceCalib.S[0].y/1000; //GE.S[1].y/1000; //pMC->S[1].y;
				det_x1 = (float)GE.DetFirstPix[1].x/1000;
				det_y1 = (float)GE.DetFirstPix[1].y/1000;
				det_x2 = (float)GE.DetLastPix[1].x/1000;
				det_y2 =  (float)GE.DetLastPix[1].y/1000;

				proj_x1 = x,
				proj_y1 = y;
				proj_x2 = (float)sliceCalib.Left[0].x/1000; //sInfo->Left[1].x/1000;
				proj_y2 =  (float)sliceCalib.Left[0].y/1000; //sInfo->Left[1].y/1000;

				sprintf(m_, "Q2:   %f, %f", x, y);
			//	MSG(m_);
							glColor3f(0.5, 0.0, 0.0);
							glBegin(GL_POLYGON);   // Source 2
							glVertex3f(x,y,0.0);
							glVertex3f(x+0.05,y,0.0);
							glVertex3f(x+0.05,y-0.05,0.0);
							glVertex3f(x,y-0.05,0.0);
							glEnd();

							glBegin(GL_LINE_STRIP);  // Detector 2
							glVertex3f(det_x1, det_y1, 0.0);
							glVertex3f(det_x2, det_y2, 0.0);
							glEnd();

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							proj_x2 = (float)sliceCalib.Right[0].x/1000; //sInfo->Right[1].x/1000;
							proj_y2 =  (float)sliceCalib.Right[0].y/1000; //sInfo->Right[1].y/1000;

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							proj_x2 = (float)sliceCalib.Ctr[0].x/1000; //sInfo->Right[0].x/1000;
							proj_y2 =  (float)sliceCalib.Ctr[0].y/1000;//sInfo->Right[0].y/1000;

							glBegin(GL_LINE_STRIP);
							glVertex3f(proj_x1, proj_y1, 0.0);
							glVertex3f(proj_x2, proj_y2, 0.0);
							glEnd();

							glBegin(GL_POLYGON);   //Center
							glVertex3f(center_Mx,center_My,0.0);
							glVertex3f(center_Mx+0.05,center_My,0.0);
							glVertex3f(center_Mx+0.05,center_My-0.05,0.0);
							glVertex3f(center_Mx,center_My-0.05,0.0);
							glEnd();

		}
		sInfo = NULL; 
	

}

GLvoid drawTomologLog(){
	GLfloat cx__, cy__, cz__;
	int i,j,r, nr_slices;
	HDC hDC;
	GLfloat radius;
	GLfloat mitte;
	float z1, z2;
	float scale;
	int alpha = LENGTH_DETAIL;  // we actually don't need to render every Slice
	static char m[256];
	fovy=12;
	float cnt;
	cnt = -0.5;

	if(myXRayLog[currentLog] != NULL){
		nr_slices = myXRayLog[currentLog]->nSlices;
		mitte = myXRayLog[currentLog]->Mess3D_Tomolog[nr_slices-1].zy/2;
		

		glPushMatrix();

			glRotatef(alfa,0.0,1.0,0.0);
			glRotatef(beta,1.0,0.0,0.0);
			glRotatef(z_angle, 0.0, 0.0, 1.0);
			
			 radius = (far_plane+near_plane)/2;

		if(defaultVModel){
			glPolygonMode(GL_FRONT, GL_FILL);
		}else if(frontVModelP) {
			glPolygonMode(GL_FRONT, GL_LINE);
		}




			for (i=0; i<nr_slices-alpha-1; i+= alpha) {
				glClear(GL_DEPTH_BUFFER_BIT);
				glPushMatrix();
				if((i > areaStart && i < areaEnd)){
						glColor3f(0.0, 1.0, 0.0);
				}else{
					glColor3f(0.5, 0.0, 0.0);
				}
					
					glTranslatef(0.0, 0.0, cnt);
					glScalef(1.0f,1.0f,1.0f);
					
					if (r<1) { r=1; }

						cx__ = (float)(myXRayLog[currentLog]->Mess3D_Tomolog[i].zx/1000);
						cy__ = (float)(myXRayLog[currentLog]->Mess3D_Tomolog[i].zy/1000);
						cz__ = (float)(myXRayLog[currentLog]->Mess3D_Tomolog[i].x/10000);
						

					LogDiskTomo(i, LOG_RADIAL_STEP, cnt, &(myXRayLog[currentLog]->Mess3D_Tomolog[i]),&(myXRayLog[currentLog]->Mess3D_Tomolog[i+(alpha)]), 
						mitte, myXRayLog[currentLog]->pGlView, XRAY);
						
					    glBegin(GL_POLYGON);
						glColor3f(1.0, 0.0, 0.0);
						glVertex3f(cx__,cy__,cz__);
						glVertex3f(cx__+0.005,cy__,cz__);
						glVertex3f(cx__+0.005,cy__-0.005,cz__);
						glVertex3f(cx__,cy__-0.005,cz__);
						glEnd();

					//	LogDisk2(i, LOG_RADIAL_STEP, cnt, &(myLogs[currentLog]->Mess3D[i]),&(myLogs[currentLog]->Mess3D[i+(alpha)]), 
					//		mitte, myLogs[currentLog]->pGlView, DISHAPE_WBARK);
					cnt -= 0.01f;
				glPopMatrix();
			}

		glPopMatrix();
		glFinish();
	//	glutSwapBuffers();
	}
}


/*
	this is the Log with Bark
  */
GLvoid drawWireframeLog(){
	GLfloat cx__, cy__, cz__, pt1x, pt1y, pt1z, pt2x, pt2y, pt2z, pt1_xo, pt1yo, pt1zo, pt2xo, pt2yo, pt2zo;;
	int i,j,r, nr_slices;
	HDC hDC;
	GLfloat radius;
	GLfloat mitte;

	float z1, z2;
	float scale;

	int alpha = LENGTH_DETAIL;  // we actually don't need to render every Slice

	fovy=12;
	float cnt;
	cnt = -0.5;


	if(myLogs[currentLog] != NULL){
		nr_slices = myLogs[currentLog]->nSlices;
		mitte = myLogs[currentLog]->Mess3D[nr_slices-1].zy/2;
		
	///	scale = nr_slices/alpha;
	//	scale /= 130000;
		glPushMatrix();

    //		glTranslatef(0.0, 0.0, -radius);
			glRotatef(alfa,0.0,1.0,0.0);
			glRotatef(beta,1.0,0.0,0.0);
			glRotatef(z_angle, 0.0, 0.0, 1.0);
			
			 radius = (far_plane+near_plane)/2;

		if(defaultVModel){
			glPolygonMode(GL_FRONT, GL_FILL);
		}else if(frontVModelP) {
			glPolygonMode(GL_FRONT, GL_LINE);
		}

		glPushMatrix();
		//	draw3DAxis();
			
		    for (i=0; i<nr_slices-alpha-1; i+= alpha) {
				glClear(GL_DEPTH_BUFFER_BIT);
				glPushMatrix();
				if(!showOnlyMark){
					if((i > areaStart && i < areaEnd)){
							glColor3f(0.7, 0.0, 0.0);
					}else{
						glColor3f(0.7, 0.7, 0.7);
					}
						
						glTranslatef(0.0, 0.0, cnt);
						glScalef(1.0f,1.0f,1.0f);
						
						if (r<1) { r=1; }
	


						cx__ = (float)(myLogs[currentLog]->Mess3D[i].zx/10000.0f);
						cy__ = (float)(myLogs[currentLog]->Mess3D[i].zy/10000.0f);
						cz__ = 0.0f; //(float)(myLogs[currentLog]->Mess3D[i].x/10000.0f);
						

						glBegin(GL_POLYGON);
						glColor3f(0.0, 0.0, 1.0);
						glVertex3f(cx__,cy__,cz__);
						glVertex3f(cx__+0.002,cy__,cz__);
						glVertex3f(cx__+0.002,cy__-0.002,cz__);
						glVertex3f(cx__,cy__-0.002,cz__);
						glEnd();

						LogDisk1(i, LOG_RADIAL_STEP, cnt, &(myLogs[currentLog]->Mess3D[i]),&(myLogs[currentLog]->Mess3D[i+(alpha)]), 
							mitte, myLogs[currentLog]->pGlView, DISHAPE_WBARK);

						//	LogDisk2(i, LOG_RADIAL_STEP, cnt, &(myLogs[currentLog]->Mess3D[i]),&(myLogs[currentLog]->Mess3D[i+(alpha)]), 
						//		mitte, myLogs[currentLog]->pGlView, DISHAPE_WBARK);
						cnt -= 0.01f;
				}else{
					glTranslatef(0.0, 0.0, cnt);
					glScalef(1.0f,1.0f,1.0f);
					
					if (r<1) { r=1; }


						cx__ = (float)(myLogs[currentLog]->Mess3D[i].zx/10000.0f);
						cy__ = (float)(myLogs[currentLog]->Mess3D[i].zy/10000.0f);
						cz__ = 0.0f; //(float)(myLogs[currentLog]->Mess3D[i].x/10000.0f);
						

						glBegin(GL_POLYGON);
						glColor3f(0.0, 0.0, 1.0);
						glVertex3f(cx__,cy__,cz__);
						glVertex3f(cx__+0.002,cy__,cz__);
						glVertex3f(cx__+0.002,cy__-0.002,cz__);
						glVertex3f(cx__,cy__-0.002,cz__);
						glEnd();

						LogDiskKr(i, LOG_RADIAL_STEP, cnt, &(myLogs[logMatch]->Mess3D[i]),&(myLogs[logMatch]->Mess3D[i+(alpha)]), mitte, myLogs[logMatch]->pGlView, DISHAPE_WBARK);


						
					if((i > areaStart && i < areaEnd)){
							glColor3f(0.7, 0.0, 0.0);
							LogDisk1(i, LOG_RADIAL_STEP, cnt, &(myLogs[logMatch]->Mess3D[i]),&(myLogs[logMatch]->Mess3D[i+(alpha)]), mitte, myLogs[logMatch]->pGlView, DISHAPE_WBARK);
					}else{
					//	glColor3f(0.0, 1.0, 0.0);
					}
					cnt -= 0.01f;
				}


	
							    pt1x = (float)(myLogs[currentLog]->Mess3D[10].zx/10000.0f);
								pt1y = (float)(myLogs[currentLog]->Mess3D[10].zy/10000.0f);
								pt1z = (float)(myLogs[currentLog]->Mess3D[10].x/10000.0f);
								pt2x = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].zx/10000.0f);
								pt2y = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].zy/10000.0f);
								pt2z = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].x/10000.0f);



						glBegin(GL_LINE_STRIP);
						glColor3f(0.0, 0.0, 1.0);
						glVertex3f(pt1x, pt1y, 0.0);
						glVertex3f(pt2x, pt2y, 0.0);
						glEnd();


								pt1x = (float)(myLogs[currentLog]->Mess3D[10].pointx[myLogs[currentLog]->degMaxKr]/10000.0f);
								pt1y = (float)(myLogs[currentLog]->Mess3D[10].pointy[myLogs[currentLog]->degMaxKr]/10000.0f);
								pt1z = (float)(myLogs[currentLog]->Mess3D[10].x/10000.0f);
								pt2x = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].pointx[myLogs[currentLog]->degMaxKr]/10000.0f);
								pt2y = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].pointy[myLogs[currentLog]->degMaxKr]/10000.0f);
								pt2z = (float)(myLogs[currentLog]->Mess3D[nr_slices-10].x/10000.0f);

						glBegin(GL_LINE_STRIP);
						glColor3f(1.0, 0.0, 0.0);
						glVertex3f(pt1x, pt1y, 0.0);
						glVertex3f(pt2x, pt2y, 0.0);
						glEnd();
	
						


				glPopMatrix();
			}







		glPopMatrix();
		glFinish();
	//	glutSwapBuffers();
	}

}

/*
	this is the Log without Bark
  */
GLvoid drawWireframeLogMatch(){
	GLfloat cx__, cy__, cz__, pt1x, pt1y, pt1z, pt2x, pt2y, pt2z;
	int i,j,r, nr_slices;
	HDC hDC;
	GLfloat radius;
	GLfloat mitte;


	int alpha = LENGTH_DETAIL;  // we actually don't need to render every Slice
	fovy=12;
	float cnt;
	cnt = -0.5;

	if(myLogs[logMatch] != NULL){
		nr_slices = myLogs[logMatch]->nSlices;
		mitte = myLogs[logMatch]->Mess3D[nr_slices-1].zy/2;
	
		glPushMatrix();
    //		glTranslatef(0.0, 0.0, -radius);
			glRotatef(alfa,0.0,1.0,0.0);
			glRotatef(beta,1.0,0.0,0.0);
			glRotatef(z_angle, 0.0, 0.0, 0.0);
			
			radius = (far_plane+near_plane)/2;

		if(defaultVModel){
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}else if(frontVModelP) {
			glPolygonMode(GL_FRONT, GL_POLYGON);
		}

		glPushMatrix();

	

			for (i=0; i<nr_slices-alpha-1; i+= alpha) {
				glClear(GL_DEPTH_BUFFER_BIT);
				glPushMatrix();
				if(!showOnlyMark){
					if((i > areaStart && i < areaEnd)){
							glColor3f(0.0, 0.0, 1.0);
					}else{
						glColor3f(0.0, 1.0, 0.0);
					}
					glTranslatef(0.0, 0.0, cnt);
					glScalef(1.0f,1.0f,1.0f);
					
					if (r<1) { r=1; }

						cx__ = (float)(myLogs[logMatch]->Mess3D[i].zx/10000.0f);
						cy__ = (float)(myLogs[logMatch]->Mess3D[i].zy/10000.0f);
						cz__ = 0.0f; //(float)(myLogs[logMatch]->Mess3D[i].x/10000.0f);;


						glBegin(GL_POLYGON);
						glColor3f(1.0, 1.0, 1.0);
						glVertex3f(cx__,cy__,cz__);
						glVertex3f(cx__+0.002,cy__,cz__);
						glVertex3f(cx__+0.002,cy__-0.002,cz__);
						glVertex3f(cx__,cy__-0.002,cz__);
						glEnd();
						LogDisk1(i, LOG_RADIAL_STEP, cnt, &(myLogs[logMatch]->Mess3D[i]),&(myLogs[logMatch]->Mess3D[i+(alpha)]), mitte, myLogs[logMatch]->pGlView, 2);

					cnt -= 0.01f;
				}else{
					glTranslatef(0.0, 0.0, cnt);
					glScalef(1.0f,1.0f,1.0f);
					
					if (r<1) { r=1; }
					if((i > areaStart && i < areaEnd)){
							glColor3f(0.0, 0.0, 1.0);
							LogDisk1(i, LOG_RADIAL_STEP, cnt, &(myLogs[logMatch]->Mess3D[i]),&(myLogs[logMatch]->Mess3D[i+(alpha)]), mitte, myLogs[logMatch]->pGlView, DISHAPE_WOBARK);
					}else{
					//	glColor3f(0.0, 1.0, 0.0);
					}
					cnt -= 0.01f;
				}

				glPopMatrix();
			}
		glPopMatrix();
		glFinish();
	//	glutSwapBuffers();
	}
}

GLvoid LogDiskTomo(int Slice_, int k, float cnt, MESS *s1,MESS *s2, float mitte,GLVIEW* pGlView, DWORD Type_){
	GLfloat p[4][2];
	GLfloat v1[3];
	GLfloat v2[3];
	GLfloat v3[3];
	GLfloat v4[3];
	GLfloat nv1[3];
	GLfloat nv2[3];
	GLfloat n[3];
	int measurePts[4];
	bool measurePoint = false;
	int i;


	measurePts[0] = DEGREE_MEASURE_D1;
	measurePts[1] = DEGREE_MEASURE_D2;
	measurePts[2] = DEGREE_MEASURE_D1+180;
	measurePts[3] = DEGREE_MEASURE_D2+180;
	
	if(s1 != NULL && s2 != NULL){ // fix


		for (i=0; i<4; i++) {
			p[0][0]=(float)s1->pointx[measurePts[i]]/10000.0f;
			p[0][1]=(float)s1->pointy[measurePts[i]]/10000.0f;

		
			if(showOnlyMark){
				if((Slice_ > areaStart && Slice_ < areaEnd)){
					if(((i > DEGREE_MEASURE) && (i < DEGREE_MEASURE+(k*2)))){
						glColor3f(0.0, 1.0, 0.0);
					}else{
							if(Type_ == DISHAPE_WBARK){
								glColor3f(0.7, 0.0, 0.0);
							}else if(Type_ == DISHAPE_WOBARK){
									glColor3f(0.0, 0.0, 1.0);
							}else if(Type_ == XRAY){
									glColor3f(0.0, 1.0, 0.0);
							}
					}

				}

			}else{

			}
					glBegin(GL_POLYGON);
					glVertex3f(p[0][0],p[0][1],0.0);
					glVertex3f(p[0][0]+0.005,p[0][1],0.0);
					glVertex3f(p[0][0]+0.005,p[0][1]-0.005,0.0);
					glVertex3f(p[0][0],p[0][1]-0.005,0.0);
					glEnd();

					if(Slice_ == areaStart){
						if(VISUALIZE_TOMO_GEOMETRY){
							drawTomoGeom(Slice_);
						}

					}

			
			

		}
	}
}


GLvoid LogDiskKr(int Slice_, int k, float cnt, MESS *s1,MESS *s2, float mitte,GLVIEW* pGlView, DWORD Type_){
	GLfloat cx__, cy__, cz__;
	GLfloat p[4][2];
	GLfloat v1[3];
	GLfloat v2[3];
	GLfloat v3[3];
	GLfloat v4[3];
	GLfloat nv1[3];
	GLfloat nv2[3];
	GLfloat n[3];
	bool measurePoint = false;
	int i;

	
	if(s1 != NULL && s2 != NULL){ // fix
		p[0][0]=(float)s1->pointx[0]/10000.0f;
		p[0][1]=(float)s1->pointy[0]/10000.0f;

		p[1][0]=(float)s2->pointx[0]/10000.0f;
		p[1][1]=(float)s2->pointy[0]/10000.0f;

		int deg = myLogs[currentLog]->degMaxKr;
		for (i=deg; i<deg+1; i+=k) {
			p[2][0]=(float)s2->pointx[deg]/10000.0f;
			p[2][1]=(float)s2->pointy[deg]/10000.0f;

			p[3][0]=(float)s1->pointx[deg]/10000.0f;
			p[3][1]=(float)s1->pointy[deg]/10000.0f;



			if(i >= myLogs[currentLog]->degMaxKr && i <= (myLogs[currentLog]->degMaxKr+LOG_RADIAL_STEP)){
			  glColor3f(1.0, 0.0, 0.0);
			}else{
							if(Type_ == DISHAPE_WBARK){
								glColor3f(0.7, 0.7, 0.7);
							}else if(Type_ == DISHAPE_WOBARK){
								glColor3f(0.0, 1.0, 0.0);
							}
			}


						cx__ = p[2][0];
						cy__ = p[2][1];
						cz__ = 0.0f; //(float)(myLogs[logMatch]->Mess3D[i].x/10000.0f);;


						glBegin(GL_POLYGON);
						glColor3f(1.0, 0.0, 0.0);
						glVertex3f(cx__,cy__,cz__);
						glVertex3f(cx__+0.002,cy__,cz__);
						glVertex3f(cx__+0.002,cy__-0.002,cz__);
						glVertex3f(cx__,cy__-0.002,cz__);
						glEnd();
			/*


					glBegin(GL_LINES);
					glVertex3f(p[0][0],p[0][1],0.0);
					glVertex3f(p[1][0],p[1][1],0.0);

					glVertex3f(p[1][0],p[1][1],0.0);
					glVertex3f(p[2][0],p[2][1],0.0);

					glVertex3f(p[2][0],p[2][1],0.0);
					glVertex3f(p[3][0], p[3][1], 0.0);

					glVertex3f(p[3][0], p[3][1], 0.0);
					glVertex3f(p[0][0],p[0][1],0.0);

					glEnd();
			
			p[0][0]=p[3][0];
			p[0][1]=p[3][1];

			p[1][0]=p[2][0];
			p[1][1]=p[2][1];
			*/
		}
	}
}




GLvoid LogDisk1(int Slice_, int k, float cnt, MESS *s1,MESS *s2, float mitte,GLVIEW* pGlView, DWORD Type_){
	GLfloat p[4][2];
	GLfloat v1[3];
	GLfloat v2[3];
	GLfloat v3[3];
	GLfloat v4[3];
	GLfloat nv1[3];
	GLfloat nv2[3];
	GLfloat n[3];
	bool measurePoint = false;
	int i;

	
	if(s1 != NULL && s2 != NULL){ // fix
		p[0][0]=(float)s1->pointx[0]/10000.0f;
		p[0][1]=(float)s1->pointy[0]/10000.0f;

		p[1][0]=(float)s2->pointx[0]/10000.0f;
		p[1][1]=(float)s2->pointy[0]/10000.0f;

		for (i=0; i<360; i+=k) {
			p[2][0]=(float)s2->pointx[i]/10000.0f;
			p[2][1]=(float)s2->pointy[i]/10000.0f;

			p[3][0]=(float)s1->pointx[i]/10000.0f;
			p[3][1]=(float)s1->pointy[i]/10000.0f;

			if(showOnlyMark){
				if((Slice_ > areaStart && Slice_ < areaEnd)){
					if(((i > DEGREE_MEASURE) && (i < DEGREE_MEASURE+(k*2)))){
						glColor3f(0.0, 1.0, 0.4);
					}else{
							if(Type_ == DISHAPE_WBARK){
								glColor3f(0.3, 0.0, 0.0);
							}else if(Type_ == DISHAPE_WOBARK){
									glColor3f(0.0, 0.0, 1.0);
							}else if(Type_ == XRAY){
									glColor3f(0.0, 0.0, 1.0);
							}
					}
				}

			}else{

			}

			if(i >= myLogs[currentLog]->degMaxKr && i <= (myLogs[currentLog]->degMaxKr+LOG_RADIAL_STEP)){
			  glColor3f(1.0, 0.0, 0.0);
			}else{
							if(Type_ == DISHAPE_WBARK){
								glColor3f(0.9, 0.9, 0.9); // was 0.7  
							}else if(Type_ == DISHAPE_WOBARK){
								glColor3f(0.0, 1.0, 0.0);
							}
			}
					glBegin(GL_LINES);
					glVertex3f(p[0][0],p[0][1],0.0);
					glVertex3f(p[1][0],p[1][1],0.0);

					glVertex3f(p[1][0],p[1][1],0.0);
					glVertex3f(p[2][0],p[2][1],0.0);

					glVertex3f(p[2][0],p[2][1],0.0);
					glVertex3f(p[3][0], p[3][1], 0.0);

					glVertex3f(p[3][0], p[3][1], 0.0);
					glVertex3f(p[0][0],p[0][1],0.0);

					glEnd();
			
			p[0][0]=p[3][0];
			p[0][1]=p[3][1];

			p[1][0]=p[2][0];
			p[1][1]=p[2][1];
		}
	}
}

GLvoid LogDisk2(int Slice_, int k, float cnt, MESS *s1,MESS *s2, float mitte,GLVIEW* pGlView, DWORD Type_){
	GLfloat p[4][2];
	GLfloat v1[3];
	GLfloat v2[3];
	GLfloat v3[3];
	GLfloat v4[3];
	GLfloat nv1[3];
	GLfloat nv2[3];
	GLfloat n[3];
	
	int fillP = 0; 
	int i;

	p[0][0]=(float)s1->pointx[0]/10000.0f;
	p[0][1]=(float)s1->pointy[0]/10000.0f;

	p[1][0]=(float)s2->pointx[0]/10000.0f;
	p[1][1]=(float)s2->pointy[0]/10000.0f;

	for (i=0; i<360; i+=k) {
		p[2][0]=(float)s2->pointx[i]/10000.0f;
		p[2][1]=(float)s2->pointy[i]/10000.0f;

		p[3][0]=(float)s1->pointx[i]/10000.0f;
		p[3][1]=(float)s1->pointy[i]/10000.0f;


		if(showOnlyMark){
			if((Slice_ > areaStart && Slice_ < areaEnd)){
				if(((i > DEGREE_MEASURE) && (i < DEGREE_MEASURE+(k*2)))){
					glColor3f(0.0, 1.0, 0.4);
				}else{
						if(Type_ == DISHAPE_WBARK){
							glColor3f(0.7, 0.0, 0.0);
						}else if(Type_ == DISHAPE_WOBARK){
								glColor3f(0.0, 0.0, 1.0);
						}

				}
			}

		}else{

		}

		glBegin(GL_POLYGON);

					v1[0]=p[1][0]-p[0][0];
					v1[1]=p[1][1]-p[0][1];
					v1[2]=0.0;

					v2[0]=p[2][0]-p[0][0];
					v2[1]=p[2][1]-p[0][1];
					v2[2]=0.0;

				normCrossProd(v1, v2, n);
				normalize(n);
				glNormal3f(n[0],n[1],n[2]);
				glVertex3f(p[0][0],p[0][1],0.0);
				glVertex3f(p[1][0],p[1][1],0.0);
				glVertex3f(p[2][0],p[2][1],0.0);
				glVertex3f(p[3][0], p[3][1], 0.0);
					
		glEnd();
			
		p[0][0]=p[3][0];
		p[0][1]=p[3][1];

		p[1][0]=p[2][0];
		p[1][1]=p[2][1];

		fillP++;
	}
}

GLvoid Disk(int k,MESS *s){
	GLfloat p[4][2];
	GLfloat v1[3];
	GLfloat v2[3];
	GLfloat n[3];
	int i;
	glBegin(GL_POLYGON);
	for (i=0; i<360; i+=k) {
		
		p[0][0]=(float)s->pointx[i]/10000.0f;
		p[0][1]=(float)s->pointy[i]/10000.0f;
	
		n[0]=0;
		n[1]=0;
		n[2]=1;
		
		glNormal3f(n[0],n[1],n[2]);
		glVertex3f(p[0][0],p[0][1],0.5);
	}

	glEnd();
}


GLvoid draw3DAxis(){
	int origin_x, origin_y;
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE);
	
	glVertex3f(0.0, 0.0, -0.5);
	glVertex3f(2.0, 2.0, -0.5);

	glVertex3f(0.0, 0.0, -0.5);
	glVertex3f(-2.0, -2.0, -0.5);
	glEnd();

}

void normalize(float v[3]){
    float d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (d != 0.0) /* avoid division by zero */
    {
		v[0] /= d;
		v[1] /= d;
		v[2] /= d;
    }
}

void normCrossProd(float v1[3], float v2[3], float out[3]){
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
    normalize(out);
}

GLvoid GetNormal(GLfloat v1[3],GLfloat v2[3],GLfloat n[3]){
	GLfloat mod;

	n[0]=v1[1]*v2[2]-v2[1]*v1[2];
	n[1]=-(v1[0]*v2[2]-v2[0]*v1[2]);
	n[2]=v1[0]*v2[1]-v2[0]*v1[1];
	mod=n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
	mod=(float)sqrt(mod);
	if (mod>0) {
		n[0]/=mod;
		n[1]/=mod;
		n[2]/=mod;
	}
	return;
}

/*
	this function does the same job as printf(..)
  */
void myPrint(float x, float y, float z, char* format, ...){
	va_list arg_list;
	char str[256];
	int i;
	va_start(arg_list, format);
	vsprintf(str, format, arg_list);
	va_end(arg_list);

	glRasterPos3f(x, y, z);
	for(i = 0; str[i] != '\0'; i++){
		glutBitmapCharacter(font_style, str[i]);
	}
}

/*
   This is where global settings are made, that is, 
   things that will not change in time 
*/
void CreateEnvironment(void){
   glEnable(GL_DEPTH_TEST);

   if (drawquality == DRAFT) {
      glShadeModel(GL_FLAT);
   }

   if (drawquality == MEDIUM) {
      glShadeModel(GL_SMOOTH);
   }

   if (drawquality == BEST) {
      glEnable(GL_LINE_SMOOTH);
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_POLYGON_SMOOTH); 
      glShadeModel(GL_SMOOTH);    
      glDisable(GL_DITHER);         /* Assume RGBA capabilities */
   }

   glLineWidth(1.0);
   glPointSize(1.0);
   glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
   glFrontFace(GL_CW);
   glDisable(GL_CULL_FACE);
   glClearColor(0.0,0.0,0.0,0.0);         /* Background colour */
   glEnable(GL_COLOR_MATERIAL);
}

/*
   This is the basic display callback routine
   It creates the geometry, lighting, and viewing position
   In this case it rotates the camera around the scene
*/
void Display(void){

	if(!LOCKED){
	   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	   glPushMatrix();
	   MakeCamera(FALSE,0,0);
	   MakeLighting();
	   MakeGeometry();
	   glPopMatrix();

	   /* glFlush(); This isn't necessary for double buffers */
	   glutSwapBuffers();
	}
}

/*
   Create the geometry
*/
void MakeGeometry(void){
   int i;
   double radius = 0.5;

	// here I draw first the board wo bark because its darker and smaller
   if(showWOBark && showWBark && showXRay){
	   //glCallList(1);
	   //glCallList(2);  // this might be more efficient, ..but the overlay is not visible
	   //glCallList(3);

	   
		  drawTomologLog();
		  drawWireframeLogMatch();
		  drawWireframeLog();
		  
	  
	
   }else{
	   if(showWOBark){
			glCallList(2);
			
	   }
		if(showWBark){
		//	glCallList(3);
			glCallList(1);
	   }
		if(showXRay){
			glCallList(3); //drawTomologVolume();
		}
   }

}


/*
   Set up the lighing environment
*/
void MakeLighting(void){
   GLfloat globalambient[] = {0.3,0.3,0.3,1.0};

   /* The specifications for 3 light sources */
   GLfloat pos0[] = {1.0,1.0,0.0,0.0};      /* w = 0 == infinite distance */
   GLfloat dif0[] = {0.8,0.8,0.8,1.0};

   GLfloat pos1[] = {5.0,-5.0,0.0,0.0};   /* Light from below */
   GLfloat dif1[] = {0.4,0.4,0.4,1.0};      /* Fainter */

   if (drawquality > DRAFT) {

      /* Set ambient globally, default ambient for light sources is 0 */
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT,globalambient);

      glLightfv(GL_LIGHT0,GL_POSITION,pos0);
      glLightfv(GL_LIGHT0,GL_DIFFUSE,dif0);

      glLightfv(GL_LIGHT1,GL_POSITION,pos1);
      glLightfv(GL_LIGHT1,GL_DIFFUSE,dif1);

      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      glEnable(GL_LIGHTING);
   }
}

/*
   Set up the camera
   Optionally creating a small viewport about 
   the mouse click point for object selection
*/
void MakeCamera(int pickmode,int x,int y){

   GLint viewport[4];

   /* Camera setup */
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (pickmode == TRUE) {
      glGetIntegerv(GL_VIEWPORT,viewport); /* Get the viewport bounds */
      gluPickMatrix(x,viewport[3]-y,3.0,3.0,viewport);
   }
   gluPerspective(Zoom,          /* Field of view */
                   1.0,          /* aspect ratio  */
                   0.1,1000.0);  /* near and far  */

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   gluLookAt(5*cos(theta*PI/180)*sin(updownrotate*PI/180),
             5*cos(updownrotate*PI/180),
             5*sin(theta*PI/180)*sin(updownrotate*PI/180), 
             0.0,0.0,0.0,                                   /* Focus    */
             0.0,1.0,0.0);                                  /* Up       */
   if (spincamera)
      theta += (cameradirection * 0.2);


   renderLog(1);
   renderLog(2);
}



/*
   Deal with plain key strokes
*/
void HandleKeyboard(unsigned char key,int x, int y){
   switch (key) {
   case 27: /* ESC */
   case 'Q':
   case 'q': exit(0); break;
   case 's':	theta += 1.0f;
   case 'S': spincamera = !spincamera; break;
   case 'b':
   case '+': Zoom -= 1.0f;  break;
   case '-': Zoom += 1.0f;  break;
   case 'r': alfa += 3.0f; break;
   case 't': beta += 3.0f; break;
   case 'z': z_angle += 3.0f; break; 
   case 'n': currentLog += 2; if(currentLog >= (MAX_LOGS3D*2)){ currentLog = 0; } renderLog(1);   renderLog(3); break;
   case 'm': logMatch = currentLog+1; if(logMatch >= MAX_LOGS3D*2){ logMatch = 0; } 
	     //  myLogs[logMatch]->pData->AlignToObject(myLogs[currentLog]->pData, (myLogs[currentLog]->nSlices/2)); 
			 myLogs[logMatch]->pData->AlignToObject(myLogs[currentLog]->pData, (myLogs[currentLog]->nSlices/2)); 
			 renderLog(2); break;
   case 'h': areaStart += 10; areaEnd = areaStart+SAMPLE_INTERVAL; renderLog(1); renderLog(3); break;
   case 'g': areaStart -=10; areaEnd = areaStart+SAMPLE_INTERVAL;   renderLog(1); renderLog(3); break;
   case 'w' : myLogs[logMatch]->pData->Rotate(10);  
	          myLogs[logMatch]->pData->AlignToObject(myLogs[currentLog]->pData, myLogs[currentLog]->nSlices/2); 
			  renderLog(2); 
			  break;
   case 'u': myXRayLog[currentLog]->sourceCalib.source2.x += 0.5f;    break;
   case 'j': myXRayLog[currentLog]->sourceCalib.source2.x -= 0.5f;    break;
   case 'i': myXRayLog[currentLog]->sourceCalib.source2.y += 0.5f;    break;
   case 'k': myXRayLog[currentLog]->sourceCalib.source2.y -= 0.5f;    break;
   case 'o': myXRayLog[currentLog]->sourceCalib.source1.x += 0.5f;    break;
   case 'l': myXRayLog[currentLog]->sourceCalib.source1.x -= 0.5f;    break;
   case 'p': myXRayLog[currentLog]->sourceCalib.source1.y += 0.5f;    break;
   case 'v': myXRayLog[currentLog]->sourceCalib.source1.y -= 0.5f;    break;
   case 'a': myXRayLog[currentLog]->calcXrayCoords(myXRayLog[currentLog], (XRAY_DIAMTYPE)0); break;
  
   case 'y':    char msg[256];   sprintf(msg, "S1_x %f  S1_y  %f        S2_x %f    S2_y  %f", myXRayLog[currentLog]->sourceCalib.source1.x, myXRayLog[currentLog]->sourceCalib.source1.y, myXRayLog[currentLog]->sourceCalib.source2.x, myXRayLog[currentLog]->sourceCalib.source2.y);
				MSG(msg);   break;
				
   }
}

/*
   Deal with special key strokes
*/
void HandleSpecialKeyboard(int key,int x, int y){
   switch (key) {
   case GLUT_KEY_LEFT:  cameradirection -= 1.0; break;
   case GLUT_KEY_RIGHT: cameradirection += 1.0;  break;
   case GLUT_KEY_UP:    updownrotate -= 1.0;  break;
   case GLUT_KEY_DOWN:  updownrotate += 1.0;  break;
   }
}

/*
   Handle mouse events
*/
void HandleMouse(int button,int state,int x,int y)
{
   int i,maxselect = 100,nhits = 0;
   GLuint selectlist[100];

   if (state == GLUT_DOWN) {
      glSelectBuffer(maxselect,selectlist);
      glRenderMode(GL_SELECT);
      glInitNames();
      glPushName(-1);

      glPushMatrix();
      MakeCamera(TRUE,x,y);
      MakeGeometry();
      glPopMatrix();
      nhits = glRenderMode(GL_RENDER);

      if (button == GLUT_LEFT_BUTTON) {

      } else if (button == GLUT_MIDDLE_BUTTON) {

      } /* Right button events are passed to menu handlers */

      if (nhits == -1)
         fprintf(stderr,"Select buffer overflow\n");

      if (nhits > 0) {
         fprintf(stderr,"\tPicked %d objects: ",nhits);
         for (i=0;i<nhits;i++)
            fprintf(stderr,"%d ",selectlist[4*i+3]);
         fprintf(stderr,"\n");
      }

   }
}

/*
   Handle the main menu
*/
void HandleMainMenu(int whichone)
{
   DWORD p;
   switch (whichone) {
   case 1: spincamera = !spincamera; break;
   case 100: exit(0); break;
   case 4: defaultVModel = true;   frontVModelP = false; break;
   case 5: frontVModelP = true;   defaultVModel = false; break;
   case 6:  findMatch(currentLog); showWBark = true; showWOBark = true;    renderLog(2);  break;
   case 7:	showWBark = true; showXRay = true; showWOBark = false;   renderLog(1); renderLog(3);   break;
   case 8: showWBark = false; showXRay = false;  showWOBark = true; renderLog(2); break;
   case 9: showWBark = true; showWOBark = true; showXRay = false; showOnlyMark = false; renderLog(100); break;
   case 10: showLogData(false);  break;   // if(!DlgLogData){    if(INFO_WIN_ID == 0){ INFO_WIN_ID = initializeInfoDlg(glutGetWindow(), 0, 250); }  }else{ DlgLogData = false;  setVisibleState(0); }
   case 13:  exportTomologVolume(); break; // showWBark = false; showWOBark = false; showXRay = true; break;
   case 20: myLogs[currentLog]->findBestMatch3D(myLogs[logMatch]);  break;
   case 30: OnShowLoadDlg(); break;
   case 40: showOnlyMark = true; renderLog(1); /*renderLog(1); renderLog(3);*/ break;
   case 41: showXRay = true;   // drawTomologVolume(); break;
   case 50: if(INFO_WIN_ID == 0){ INFO_WIN_ID = initialize3DInfoDlg(glutGetWindow(), 0, 250); } break;
   case 66: VisualizeKruemmung(); myLogs[currentLog]->lineKr.krAvail = true; break;
   	

   }
}


/*
	Here we just repeat the steps from the menu;
*/
void perform3DModelAnalysis(){
	
	if(!TRUST_SLX_MATCHING){
		findMatch(currentLog);
	}else{
		logMatch = currentLog+1; // in this case we just take te succsessive entry
	}
	
	myLogs[currentLog]->setLock(true);
	myLogs[currentLog]->findBestMatch3D(myLogs[logMatch]);
	
}

/*
   Handle the ball speed sub menu
*/
void HandleSpeedMenu(int whichone)
{
   switch (whichone) {
   case 1: break;
   case 2:    break;
   case 3:   break;
   }
}

/*
   How to handle visibility
*/
void HandleVisibility(int visible)
{
   if (visible == GLUT_VISIBLE)
      glutIdleFunc(HandleIdle);
   else
      glutIdleFunc(NULL);
}

/*
   What to do on an idle event
*/
void HandleIdle(void)
{

   glutPostRedisplay();
}

/*
   Draw text in the x-y plane
   The x,y,z coordinate is the bottom left corner (looking down -ve z axis)
*/
void DrawTextXY(double x,double y,double z,double scale,char *s)
{
   int i;

   glPushMatrix();
   glTranslatef(x,y,z);
   glScalef(scale,scale,scale);
   for (i=0;i<strlen(s);i++)
      glutStrokeCharacter(GLUT_STROKE_ROMAN,s[i]);
   glPopMatrix();
}


/*
   Display the program usage information
*/
void GiveUsage(char *cmd)
{
   fprintf(stderr,"Tomolog3D");
   exit(-1);
}


long WINAPI MatchFunc(Model3D *obj){
				obj->threadRun = true;
				int currdist, len, deg, deg2, mindist__, min_distDiam__, mindist__2, lenMM;
				double d_currdist;
				int min_dist = 1000000000;
				int min_diamDiff = 1000000000;

				int curr_diamDiff = 0;

				int min_distTmp = min_dist;
				int min_dist2 = min_dist;
				char msg_[36];

				int m25_len, m50_len, m75_len;
				
				deg = 0;
				len = obj->nSlices/2;  // here we consider slice-weise
				lenMM = obj->pData->GetLogLength()/2;

				m50_len = len;
				m25_len = len/2;
				m75_len = len+m25_len;
				int loopMeas[3];
				int loopDeg[3];

				// we split the matching accross 3 parts
				loopMeas[0] = m50_len; 
				loopMeas[1] = m25_len;
				loopMeas[2] = m75_len;

				mindist__ = min_distDiam__ = 0;

				if(obj != NULL){
					for(int c = 0; c < 1; c++){
						while(obj->threadRun){
							currdist = (int)obj->CalcDist(&(obj->pModelMatch->Mess3D[loopMeas[c]]), &(obj->Mess3D[loopMeas[c]]));
						//	d_currdist = CalcDist(&(obj->pModelMatch->Mess3D[loopMeas[c]]), &(obj->Mess3D[loopMeas[c]]));
							min_dist = min(min_dist, currdist);
							
							//we match also against tomolog ST2
						//	curr_diamDiff = abs((getLogDiam_WBark(lenMM, DEGREE_MEASURE_D1)-getLogDiam_WBarkTOMOLOG(lenMM, DEGREE_MEASURE_D1)));  // we use ST2 alg
						//	min_diamDiff = min(curr_diamDiff, min_diamDiff);

			
							if(min_dist == currdist){
								mindist__ = deg;
							}

							
							if(min_diamDiff == curr_diamDiff){
								min_distDiam__ = deg;
							}
							

							obj->pModelMatch->pData->Rotate(ROTATION_FAKTOR);
							obj->pModelMatch->pData->AlignToObject(obj->pData, len);

							if(!COLLECT_DATA){
								renderLog(2); // we render the object again, after rotation
								Sleep(10);
							}
							deg += ROTATION_FAKTOR;
							deg2 += ROTATION_FAKTOR;

							if(deg > 360){
								obj->threadRun = FALSE;
								break;
							}
						
						}
						deg = 0;
						loopDeg[c] = mindist__;
						obj->threadRun = TRUE;  // we repeat it for different parts of the object
 					}

					sprintf(msg_, "degDiShape1: %d   degTOMO: %d ", mindist__, min_distDiam__);
					MSG(msg_);




					obj->pModelMatch->pData->Rotate(mindist__); 
					obj->pModelMatch->pData->AlignToObject(obj->pData, len);

					renderLog(2);

					obj->logInfo_.RecFactor = min_dist/100000;
					obj->logInfo_.RecFactor /= 1000;
					  

				}
				obj->setLock(false);   // this is just to keep a syncronization between the loading module
				return(TRUE);
}

void OnShowLoadDlg() {
   CString TargetPath;
   CString LastBrowsedPath;
   const XRAY_SLICEINFO *pInfo; 
   FILE *pF;
   int obj_cnt = 0;
   int log_cnt = 1;
   int DiShape_offset = 0;
   char msg_[256];
   int t__; 

   	pXR = new CXRayData(); 

   if( LastBrowsedPath.IsEmpty() ){
        LoadObjDlg.m_sBrowsePath = TargetPath;
   }else{
        LoadObjDlg.m_sBrowsePath = LastBrowsedPath;
   }

    //create a new object to hold the result:
    CStringList* pStringList = new CStringList;

    LoadObjDlg.SetResultList( pStringList );

    if( LoadObjDlg.DoModal() == IDOK )
    {
        //replace the old list:
        delete m_pFileToLoadList;
        m_pFileToLoadList = pStringList;

        m_CurrentPosition = m_pFileToLoadList->GetHeadPosition(); 

		for(m_CurrentPosition = m_pFileToLoadList->GetHeadPosition(); m_CurrentPosition != NULL; ){
			mif_file =  m_pFileToLoadList->GetAt( m_CurrentPosition );
			
			if(!COLLECT_DATA && (obj_cnt < MAX_LOGS3D)){
				if(!LOAD_ONLY_XRAY){ 
					if( !LoadDataObject( m_pFileToLoadList->GetAt( m_CurrentPosition ), &LoadMifObj[obj_cnt]) ){
							MessageBox(NULL, "Unable to load .mif file!!", "ERROR", MB_OK);
       
					 }else{
						LOCKED = true; // we lock the engine in order to avoid access violations
						p3D1 = &LoadMifObj[obj_cnt].DataShape1; 
						p3D2 = &LoadMifObj[obj_cnt].DataShape2;
						pXRay = &LoadMifObj[obj_cnt].DataXRay;
				
						if((p3D1 && p3D2 && pXRay) != NULL){
							// we populate up to MAX_LOGS3D*2
							createObj3D(p3D1->get3DSlices(), p3D1->getNr3DSlices(), NULL, p3D1, DiShape_offset);
							createObj3D(p3D2->get3DSlices(), p3D2->getNr3DSlices(), NULL, p3D2, DiShape_offset+1);   DiShape_offset += 2;
							createTomologObject3D((CXRayData*)pXRay, p3D1->getNr3DSlices(), 0);

							setCurrentObj(0);

						//	exportTomologVolume();

							t__ = p3D1->IsTipForward();
							MessageBox(NULL, "Obj loaded...", "Info", MB_OK);
							sprintf(msg_, "Zopf vorne: %d", t__);
							MSG(msg_);
							obj_cnt++;
							
						}
						LOCKED = false;
					 }
				}else{   

							static char str[256];
							static char line[256];
							static char msg[256];
							const char *str1, *str2;
							CString file_name;
							FILE *pF;
							CTomologUtil *pUtil;

							float stdDiam, D1, D2, diam_woB_D1, diam_woB_D2;
							int size, validSliceInd, bValid, iD1, iD2, Sx1, Dx1, Sx2, Dx2, iMuMax_1, iMuD_1, iMuD2_1, iST_1, iST2_1, iMuMax_2, iMuD_2, iMuD2_2, iST_2, iST2_2, AVG_ALG_1, AVG_ALG_2 = 0;

							if(pXR){
				
								pXR->LoadFromTomologMif(mif_file);
								pXR->BarkAnalysis();

								if(EXPORT_DENSITY_PROFILES){
									pF = fopen(ANALYSIS_DATA_FILE_DENSITY, "a");
								}else{
									pF = fopen(ANALYSIS_DATA_FILE_EICH, "a");
								}

								if(pF){

									sprintf(msg_, "%s;\n", mif_file);
									size = strlen(msg_);
									fwrite(msg_, 1, size, pF);
									fflush(pF);


									if(!EXPORT_DENSITY_PROFILES){
										for(int i = 65; i < 90; i++){
											pInfo = pXR->GetSliceInfoByNumber(i);
											if(pInfo){
												if(pInfo->bValid){
													if(!CALC_XRAX_BARK_DETECT){
														iD1 = (int)(pInfo->D1*100);
														iD2 = (int)(pInfo->D2*100);
														Sx1 = (int)(pInfo->Sx[0]);
														Dx1 = (int)(pInfo->Dx[0]);
														Sx2 = (int)(pInfo->Sx[1]);
														Dx2 = (int)(pInfo->Dx[1]);
														sprintf(msg_, "%d;%d;%d;%d;%d;%d;%d;\n",i, Sx1, Dx1, Sx2, Dx2, iD1, iD2);
													}else{
														iD1 = (int)(pInfo->D1*100);
														iD2 = (int)(pInfo->D2*100);

														iMuMax_1 = (int)((pInfo->RSxMuMax[0]+pInfo->RDxMuMax[0])*100);
														iMuD_1 = (int)((pInfo->RSxMuD[0]+pInfo->RDxMuD[0])*100);
														iMuD2_1 = (int)((pInfo->RSxMuD2[0]+pInfo->RDxMuD2[0])*100);
														iST_1 = (int)((pInfo->RSxST[0]+pInfo->RDxST[0])*100);
														iST2_1 = (int)((pInfo->RSxST2[0]+pInfo->RDxST2[0])*100);

														AVG_ALG_1 = (iMuMax_1+iMuD_1+iMuD2_1+iST_1+iST2_1)/5;


														iMuMax_2 = (int)((pInfo->RSxMuMax[1]+pInfo->RDxMuMax[1])*100);
														iMuD_2 = (int)((pInfo->RSxMuD[1]+pInfo->RDxMuD[1])*100);
														iMuD2_2 = (int)((pInfo->RSxMuD2[1]+pInfo->RDxMuD2[1])*100);
														iST_2 = (int)((pInfo->RSxST[1]+pInfo->RDxST[1])*100);
														iST2_2 = (int)((pInfo->RSxST2[1]+pInfo->RDxST2[1])*100);

														AVG_ALG_2 = (iMuMax_2+iMuD_2+iMuD2_2+iST_2+iST2_2)/5;

														sprintf(msg_, "%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;\n",i, iD1, iD2, iMuMax_1, iMuMax_2, iMuD_1, iMuD_2, iMuD2_1, iMuD2_2, iST_1, iST_2, iST2_1, iST2_2, AVG_ALG_1, AVG_ALG_2);
													}

													if(!isMeasureOverrun(pInfo)){  // we just write acceptable results to file
														size = strlen(msg_);
														fwrite(msg_, 1, size, pF);
														fflush(pF);
													}




												}
											}
										}
									//	createTomologObject3D((CXRayData*)pXR, 150, 0);  // with 150 rows (slices)
										setCurrentObj(0);

									}else{
										const XRAY_HEADER *pHead = NULL;
										int radT = 0;
										float rad = 0;
										float radDensity = 0;
										float rsX, rdX = 0;
										const float *pProj;
										int sxIndex, dxIndex = 0;  // indexes of the projection
										
										pHead = pXR->ExtractHeaderInfo();
										DWORD   TotalLines    = pHead->NumOfRows;
										int     ByteLineWidth = pHead->RowByteSize;
										int     NumPixel      = pHead->NumPixels;
		
										for(int q = 0; q < 2; q++){  // Source iteration
											pProj = (pXR)->GetProjection(q);
											for(int i = 10; i < 190 ; i++){
												pInfo = pXR->GetSliceInfoByNumber(i);
												if(pInfo){
													if(pInfo->bValid){
															rsX = pInfo->RSx[q];
															rdX = pInfo->RDx[q];

															sxIndex = pInfo->Sx[q];
															dxIndex = pInfo->Dx[q];

															sprintf(msg_, "%d;%d;%d;%d;%d;%f;\n", log_cnt, i, q, radT, rad, radDensity);

														size = strlen(msg_);
														fwrite(msg_, 1, size, pF);
														fflush(pF);


													}
												}
											}
										}


										log_cnt++;


									}


									sprintf(msg_, "\n\n\n");
									size = strlen(msg_);
									fwrite(msg_, 1, size, pF);
									fflush(pF);
									fclose(pF);

									
								//	createTomologObject3D((CXRayData*)pXR);

								//	setCurrentObj(0);

									/*
									sprintf(msg_, "Create 3D Model ACIS-Alg?");
									MSG(msg_);

									createTomologObject3D((CXRayData*)pXR);

									setCurrentObj(0);

									exportTomologVolume();
									*/


								}else{
										MSG("log3DEich.csv not Found!!");
					
								}
					
								
									/*
								pUtil = new CTomologUtil();
								pUtil->setXRayData(pXR);		
								pUtil->DoModal();
								*/
								
								
								
							}
			
				}

			}else{// in this case we fill a fileList and than start a Thread for the analysis
				fileList[obj_cnt] = mif_file;  obj_cnt++;
			}

			 m_pFileToLoadList->GetNext(m_CurrentPosition);
        }
		
    }
    else
    {
        //free the memory we allocated, as we don't use it.
        delete pStringList;
    }	
	if(COLLECT_DATA){
		pModelAnalysis = new Model3DAnalysis(obj_cnt);
		
	}

}


void OnLoadMifFiles(){
	int cnt = 0;
	std::string s_tmp;
	int beg = 1;

   oslink::directory dir(MIF_FILEDIR);
   while (dir){
	   s_tmp = dir.next();
	   if(!beg){
			fileList[cnt].Format("%s\\%s", MIF_FILEDIR, s_tmp.c_str()); 
			showfileList[cnt].Format("%s", s_tmp.c_str());   cnt++;
	   }
	   beg = 0;
	}

	if(COLLECT_DATA){
		pModelAnalysis = new Model3DAnalysis(cnt);
	}
}

void Model3DAnalysis::run(){
   char msg_[256];
   CXRayData *pXR;
   CMirData *pMir;
   CPhotoData *pPhoto;

   /*  in pMir->m_LogData.stammdaten.
	int idMessung1;
	int	idMessung2;
	int idXray;
    int idMir;
	int idScreenlog;

	WORD fErr;
	WORD flags;
	WORD auxflags;

	char ho [KUERZ+1] ;			// Holzart
	char qk [KUERZ+1] ; //übernahme   // Qualität
	char qs [KUERZ+1] ;  //sortierung	// Qualität Sortiment
	
	int la;			// Länge
	int zo;			// Zopf
	int rmi;		// Durchmesser
	int abh;			// Abholzigkeit
	int kr;				// Krümmung

	char klasse[KUERZ+1];       //
	char klassifizierung[13];   
	char sortierung[13];

	char mkla[KLA_LENGTH+1];   //klasse mittendurchm.
	
	int kla;        //kapp länge
	int fm;
	int box;
	int sort;
	int ko;
	int rt;			// RindenTabelle

  */
	int obj_cnt = 0;
   int log_cnt = 1;//1486
	while(bThreadRunning){
		for(int i = 0; i < nr_files; i++){
				if( !LoadDataObject(fileList[i], &LoadMifObj[obj_cnt]) ){
						MessageBox(NULL, "Unable to load .mif file!!", "ERROR", MB_OK);
						if(LOAD_DISHAPE3D_MODELS){
							goto BEGINN;   // not elegant, but necessary
						}
				 }else{
					if(!KR_ANALYSIS){
						LOCKED = true; // we lock the engine in order to avoid access violations
						p3D1 = &LoadMifObj[obj_cnt].DataShape1; 
						p3D2 = &LoadMifObj[obj_cnt].DataShape2;
						pXRay = &LoadMifObj[obj_cnt].DataXRay;
						pMir = &LoadMifObj[obj_cnt].DataMir;
				
						if((p3D1 && p3D2 && pXRay && pMir) != NULL){
							createObj3D(p3D1->get3DSlices(), p3D1->getNr3DSlices(), NULL, p3D1, 0);
							createObj3D(p3D2->get3DSlices(), p3D2->getNr3DSlices(), NULL, p3D2, 1);
							createTomologObject3D((CXRayData*)pXRay, p3D1->getNr3DSlices(), 0);

							setCurrentObj(0);  // set currentLog = 0;


							sprintf(myLogs[currentLog]->mirData.ho, "%s", pMir->m_LogData.stammdaten.ho);
							sprintf(myLogs[currentLog]->mirData.qk, "%s", pMir->m_LogData.stammdaten.qk);
							sprintf(myLogs[currentLog]->mirData.qs, "%s", pMir->m_LogData.stammdaten.qs);
							myLogs[currentLog]->mirData.la = pMir->m_LogData.stammdaten.la;
							myLogs[currentLog]->mirData.la = pMir->m_LogData.stammdaten.zo;
							myLogs[currentLog]->mirData.rmi = pMir->m_LogData.stammdaten.rmi;
							myLogs[currentLog]->mirData.abh = pMir->m_LogData.stammdaten.abh;
							myLogs[currentLog]->mirData.kr = pMir->m_LogData.stammdaten.kr;


							myLogs[currentLog]->logInfo_.mif_file = showfileList[i];  // we store the file name
							myLogs[currentLog]->logInfo_.ID = log_cnt; log_cnt++;

							perform3DModelAnalysis();	// in here a matching thread is started, therefore we have to wait for the lock
							while(myLogs[currentLog]->getLock()){  // here we wait for the obj to complete the operation
								Sleep(2000);
							}
							showLogData(true); // we print the results to File

							obj_cnt = 0;		// in this case we always fill the same obj;) 
							 


							//we free our resources
							/*
							myLogs[currentLog]->~Model3D();
							myLogs[logMatch]->~Model3D();
							myXRayLog[currentLog]->~Model3D();

							free(myXRayLog[currentLog]);
							free(myLogs[currentLog]);
							free(myLogs[logMatch]);
							*/
							
						}else{
							log_cnt++;
						}

					}else{  // for curvature Analysis (Krümmung) KR_ANALYSIS == 1
BEGINN:
						if(LOAD_DISHAPE3D_MODELS){  // in this case we just load the DiShape struct
							memset(&myL, 0, sizeof(myL));
						//	ReadMessung((LPTSTR)fileList[i].GetBuffer(200), &myL);
							createObj3D(NULL, 0, NULL, NULL, 0); // we just have to avoid NULLS for exportLogDataKr(bool);
							createObj3D(NULL, 0, NULL, NULL, 1);
							myLogs[currentLog]->logInfo_.mif_file = showfileList[i];
						//	exportLogDataKr2(true);


						}else{
							LOCKED = true; // we lock the engine in order to avoid access violations
							p3D1 = &LoadMifObj[obj_cnt].DataShape1; 
							p3D2 = &LoadMifObj[obj_cnt].DataShape2;
							pMir = &LoadMifObj[obj_cnt].DataMir;
							pPhoto = &LoadMifObj[obj_cnt].DataSL;
							if((p3D1 && p3D2 && pMir) != NULL){
								createObj3D(p3D1->get3DSlices(), p3D1->getNr3DSlices(), NULL, p3D1, 0);
								createObj3D(p3D2->get3DSlices(), p3D2->getNr3DSlices(), NULL, p3D2, 1);
								

								if(!SIMULATE_WZLR_SIGNAL){
									if(pMir->m_LogData.stammdaten.reserved[0] == 1){  //reduziert
											setCurrentObj(0);  // we take 3D1   with Bark
									}else{
											setCurrentObj(1); // we take 3D2   without Bark
									}
								}else{  // we check if the log would have been reduces
									int diam1, diam2, diam3, diam11, diam22, diam33, len;
									len = p3D1->GetLogLength()*10;  //mm
									if(p3D1->IsTipForward()){
										  diam1 = p3D1->GetDiameterAt(len-50, 100, false);
										  diam2 = p3D1->GetDiameterAt(len-200, 100, false);
										  diam3 = p3D1->GetDiameterAt(len-500, 100, false);

										  diam11 = p3D1->GetDiameterAt(len-50, 10, false);
										  diam22 = p3D1->GetDiameterAt(len-200, 10, false);
										  diam33 = p3D1->GetDiameterAt(len-500, 10, false);
										
									}else{
										  diam1 = p3D1->GetDiameterAt(50, 100, false);
										  diam2 = p3D1->GetDiameterAt(200, 100, false);
										  diam3 = p3D1->GetDiameterAt(500, 100, false);

										  diam11 = p3D1->GetDiameterAt(50, 10, false);
										  diam22 = p3D1->GetDiameterAt(200,10, false);
										  diam33 = p3D1->GetDiameterAt(500, 10, false);
									}

									if(abs(diam1-diam3) >= MAX_WZLR_DIFF){
											setCurrentObj(0); // we choose 3D1
									}else{
											setCurrentObj(1);	// we take 3D2
									}

									if(abs(diam11-diam33) >= MAX_WZLR_DIFF){
											setCurrentObj(0); // we choose 3D1
									}else{
											setCurrentObj(1);	// we take 3D2
									}

								}

								if(ALWAYS_TAKE_3D1_MEASURE){
									setCurrentObj(0);
								}


								myLogs[currentLog]->pPhoto = pPhoto;  //  pointer to the Photo-Data

								sprintf(myLogs[currentLog]->mirData.ho, "%s", pMir->m_LogData.stammdaten.ho);
								sprintf(myLogs[currentLog]->mirData.qk, "%s", pMir->m_LogData.stammdaten.qk);
								sprintf(myLogs[currentLog]->mirData.qs, "%s", pMir->m_LogData.stammdaten.qs);
								myLogs[currentLog]->mirData.la = pMir->m_LogData.stammdaten.la;
								myLogs[currentLog]->mirData.zo = pMir->m_LogData.stammdaten.zo;
								myLogs[currentLog]->mirData.rmi = pMir->m_LogData.stammdaten.rmi;
								myLogs[currentLog]->mirData.abh = pMir->m_LogData.stammdaten.abh;
								myLogs[currentLog]->mirData.kr = pMir->m_LogData.stammdaten.kr;

								myLogs[currentLog]->logInfo_.mif_file = showfileList[i];  // we store the file name
								myLogs[currentLog]->logInfo_.ID = log_cnt; log_cnt++;
								

								exportLogDataKr(true); // here we don't need synchronization as above, because here we don't have a matching thread


								if(COLLECT_DATA){
								// we destroy the created obj
									myLogs[0]->~Model3D();
									myLogs[1]->~Model3D();

									free(myLogs[0]);
									free(myLogs[1]);
								}else{   // we turn the log with the kr downwards




								}

								obj_cnt = 0;		// in this case we always fill the same obj;) 
							}else{
								log_cnt++;
							}

						}
					}

				 }
		}
		bThreadRunning = false;
		sprintf(msg_, "Model3D Analysis finished");
		MSG(msg_);
		this->~Model3DAnalysis();

	}
}

typedef struct {
	double x, y, z;
}XYZ;

XYZ interpolation(XYZ p1,XYZ p2,XYZ p3,XYZ p4,double mu)
{
   double mu2=mu*mu,mu3=mu*mu*mu;
   XYZ p;

   p.x = 0.5 * ((-p1.x + 3*p2.x -3*p3.x + p4.x)*mu3 + (2*p1.x -5*p2.x + 4*p3.x - p4.x)*mu2 + (-p1.x+p3.x)*mu + 2*p2.x);
   p.y = 0.5 * ((-p1.y + 3*p2.y -3*p3.y + p4.y)*mu3 + (2*p1.y -5*p2.y + 4*p3.y - p4.y)*mu2 + (-p1.y+p3.y)*mu + 2*p2.y);
   p.z = 0.5 * ((-p1.z + 3*p2.z -3*p3.z + p4.z)*mu3 + (2*p1.z -5*p2.z + 4*p3.z - p4.z)*mu2 + (-p1.z+p3.z)*mu + 2*p2.z);

   return p;
}



void InterpolateValues(int *x__, int *y__, int size__, int *res_x, int *res_y){
	static int MAX_RESOLUTION = 1000;
	int i, cnt;

	XYZ		inp[ 100 ];
	int		min, max;

	cnt = 0;
	min = MAX_RESOLUTION;
	max = 0;


	for ( i = 0; i < size__; i++ ) {
			inp[ i ].x = (float)x__[i];   *x__++;
			inp[ i ].y = (float)y__[i];    *y__++;
			inp[ i ].z = 0;
			cnt++;

	}


		int res;
		int j, k, j2, n;
		double mu, dmu;
		int p1,p2,p3,p4;
		XYZ o;

		n = cnt;
		k = min;
		j2 = 0;
		mu = 0;
		res = MAX_RESOLUTION;
		dmu = 1.0*n/res;

		for ( i = 0; i < res; i++ ) {
			j = i*n/res;
			if ( j>j2) mu = 0;	
			p1 = max(0,j-1);
			p2 = j;
			p3 = min(n-1,j+1);
			p4 = min(n-1,j+2);
			o = interpolation( inp[p1],inp[p2],inp[p3],inp[p4], mu );
			*res_x = (int)o.x;   *res_x++;
			*res_y = (int)o.y;	*res_y++;
			mu+=dmu;
			mu=min(1.0,mu);
			j2=j;
		}
		

}


/*
	we could solve this by taking the derivative of oder -n, but we simplify it to a line segment alg
  */
int hasKrEinseitig(LOG *plog__, int z1, int z2){
	/*
	float d;
	float *pCurv;

	int in_x[2];
	int in_y[2];
	int res_x[1000];
	int res_y[1000];
	int f_, x_, y_, point_1_x, point_2_x, point_1_y, point_2_y, c1x, c1y, c2x, c2y, c3x, c3y, logCtr;

	point_1_x = plog__->slices[z1].centerx;
	point_1_y = plog__->slices[z1].centery;

	point_2_x = plog__->slices[z2].centerx;
	point_2_y = plog__->slices[z2].centery;

	in_x[0] = point_1_x;
	in_x[1] = point_2_x;
	in_y[0] = point_1_y;
	in_y[1] = point_2_y;




	InterpolateValues(&in_x[0], &in_y[0], 2, &res_x[0], &res_y[0]);
	int cnt = 0;

	c1x = plog__->slices[z1+1].centerx;
	c2x = plog__->slices[z2-1].centerx;

	c1y = plog__->slices[z1+1].centery;
	c2y = plog__->slices[z2-1].centery;

	int dir_1, dir_2;
	bool einseitig = 0;

	int horr, vert = 0;
	if(point_1_x > point_2_x || point_1_x < point_2_x){
			vert = 1;
	}else if(point_1_y > point_2_y || point_1_y < point_2_y){
			horr = 1;
	}

	int tmp_dir = 0;
	dir_1 = 0;
	for(int i = z1+1; i < z2-1; i++){
			c1x = plog__->slices[i].centerx;
			c1y = plog__->slices[i].centery;


			
			dir_1 = c1x-res_x[cnt];
			dir_2 = c1y-res_y[cnt];

			if(horr){
				if(c1x > res_x[cnt]){  // horrizontal
					dir_1 = 1;
					tmp_dir = dir_1;
				}

				if(c1x > res_x[cnt]){  // horrizontal
					dir_1 = 1;
					tmp_dir = dir_1;
				}
			}

			if(vert){
			   	if(c1y > res_y[cnt]){  // horrizontal
					dir_1 = 1;
					tmp_dir = dir_1;
				}
			}

			cnt++;
	}

  */
	return 1;

}


void smoothLog(Model3D *pModel){
	int median_x = 0;
	int median_y = 0;
	int windowSize = SMOOTH_WIN_SIZE;
	int nr = 0;
	int min_x[SMOOTH_WIN_SIZE];
	int min_y[SMOOTH_WIN_SIZE];

	int off_s = 10;
	int off_e = 0;
	
	for(int deg = 0; deg <= 360; deg++){
	//	if(deg % 10 == 0){   		renderLog(currentLog+1);  }
		nr = 0;
		for(int i = 10; i < (pModel->nSlices-windowSize)-10; i++){
			if((i != 10) && ((i % windowSize) == 0)){
				for(int x = off_s; x <= windowSize; x++){
					if(x < pModel->nSlices){
						// we find the min
							qsort(min_x, windowSize, sizeof(int), compareMin);
							qsort(min_y, windowSize, sizeof(int), compareMin);

						pModel->Mess3D[x].pointx[deg] = (int)min_x[0]; //(median_x/windowSize);
						pModel->Mess3D[x].pointy[deg] = (int)min_y[0]; //(median_y/windowSize);		 
					}
				}
				off_s += windowSize;
				median_x = 0;
				median_y = 0;

			}else{
				median_x += pModel->Mess3D[i].pointx[deg];
				median_y += pModel->Mess3D[i].pointy[deg];

				min_x[i%windowSize] = pModel->Mess3D[i].pointx[deg];
				min_y[i%windowSize] = pModel->Mess3D[i].pointy[deg];
			}


		}
	}

}

int getOvalitaet(int degKr, int pos){
	/*
	int d1, d2, op;  // d1 is the diameter in degMaxKr
	degkr = (360-degKr)-90;
	op = (degKr >= 90) ? (degKr-90) : (degKr+90);
	d1 = myLogs[currentLog]->pData->GetDiameterAt(pos, degKr, bClambDiShapeDiam);
	d2 = myLogs[currentLog]->pData->GetDiameterAt(pos, op, bClambDiShapeDiam);
	return (d1-d2);
	*/
	return 0;
}

void exportLogDataKr(bool writeToFile){
	static int offLenZopf[100];
	static int offLenKopf[100];
	static int values[25][25];
	static int valuesMantel[25][25];
	static char img_path[256];
	double winkelKr = 0;
	int einseitigeKr;
	int cntWA = 0;
	int krTmp, krTmp2, krTmp3, krTmp4;
	CString S1 = "S1";
	CString S2 = "S2";

	int off_start = 0;
	int off_end = 0;
	memset(&dataSetKr, 0, sizeof(KRData));
	memset(&offLenZopf, 0, sizeof(offLenZopf));
	memset(&offLenKopf, 0, sizeof(offLenKopf));
	memset(&myL, 0, sizeof(LOG));
	memset(&values, 0, sizeof(values));
	memset(&valuesMantel, 0, sizeof(valuesMantel));


	int degMaxKr = 0;
	// len offsets in % from Zopf
	offLenZopf[0] = 3;
	offLenZopf[1] = 5;
	offLenZopf[2] = 8;
	offLenZopf[3] = 100; // this is indeed 10 cm  -> 
				
	// len offsets in % from Kopf
	offLenKopf[0] = 6;
	offLenKopf[1] = 8;
	offLenKopf[2] = 10;
	offLenKopf[3] = 12;
	offLenKopf[4] = 16;
	offLenKopf[5] = 20;

	if(myLogs[currentLog] != NULL && myLogs[logMatch] != NULL){
		sprintf(dataSetKr.mif_file, "%s", myLogs[currentLog]->logInfo_.mif_file.GetBuffer(256));
		
		myL.la = myLogs[currentLog]->getLength()*10;  // we procceed in mm
		myL.nrslices = myLogs[currentLog]->nSlices;
		myL.lage =   (myLogs[0]->pData)->IsTipForward(); // Zopf priori   here we always take 3D1 as measure  myLogs[0]

		int diamRMI = myLogs[currentLog]->mirData.rmi; 

	//	smoothLog(myLogs[currentLog]);   // we smooth the surface
		for(int i = 0;i < myLogs[currentLog]->nSlices; i++){
			myL.slices[i].centerx = myLogs[currentLog]->Mess3D[i].zx;
			myL.slices[i].centery = myLogs[currentLog]->Mess3D[i].zy;
			myL.slices[i].z = myLogs[currentLog]->Mess3D[i].x;

			for(int j =0; j < 360; j++){
				myL.slices[i].fpoint[j].x = myLogs[currentLog]->Mess3D[i].pointx[j];
				myL.slices[i].fpoint[j].y = myLogs[currentLog]->Mess3D[i].pointy[j];
			}
			/*
			int sx = 0;
			int sy = 0;
			for(int j =0; j < 360; j++){
				myL.slices[i].fpoint[j].x = myLogs[currentLog]->Mess3D[i].pointx[j];
				myL.slices[i].fpoint[j].y = myLogs[currentLog]->Mess3D[i].pointy[j];
				sx+=myL.slices[i].fpoint[j].x;
				sy+=myL.slices[i].fpoint[j].y;
			}
			myL.slices[i].centerx=sx/360;			// we also need to adjust the central point
			myL.slices[i].centery=sy/360;
			*/
		}



		dataSetKr.KR_EINSEITIG = hasKrEinseitig(&myL, 10, myL.nrslices-10);

		// we reset the offsets
		off_end = 0;
		off_start = 0;
		
		myL.einseitigeKr_90  = 0;
		myL.einseitigeKr_180 = 0;
		

		// this algorithm is used to detect a "Einseitiger-Wurzelalnauf"		
		int steigungB, steigungD, begX, begY, endX, endY, ptA, ptB = 0;
		int invertedIndex_1, invertedIndex_2 = 0;
		if(myL.lage){ // Zopf at beginning
			invertedIndex_1 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF;
			invertedIndex_2 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_B;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery
				;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungB = ptA/ptB;
			}else{
				steigungB = ptA/1;
			}

			invertedIndex_1 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_C;
			invertedIndex_2 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_D;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungD = ptA/ptB;
			}else{
				steigungD = ptA/1;
			}

			if(steigungB > steigungD+DIFF_OFFSET_WURZEL_ANLAUF){
				dataSetKr.WA = 1;
			}else{
				dataSetKr.WA = 0;
			}
		}else{	// Zopf at end
			invertedIndex_1 = OFFSET_WURZEL_ANLAUF;
			invertedIndex_2 = OFFSET_WURZEL_ANLAUF_B;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungB = ptA/ptB;
			}else{
				steigungB = ptA/1;
			}


			invertedIndex_1 = OFFSET_WURZEL_ANLAUF_C;
			invertedIndex_2 = OFFSET_WURZEL_ANLAUF_D;

			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungD = ptA/ptB;
			}else{
				steigungD = ptA/1;
			}

			if(steigungB > steigungD+DIFF_OFFSET_WURZEL_ANLAUF){
				dataSetKr.WA = 1;
			}else{
				dataSetKr.WA = 0;
			}

		}

		// we store both WA-factors in the Dataset
		dataSetKr.wurzelAn.stB = steigungB;
		dataSetKr.wurzelAn.stD = steigungD;
		dataSetKr.wurzelAn.stDiff = steigungB-steigungD;  // if negative -> steigungD > steigungB
		


		for(int c = 0; c < 4; c++){
			for(int j = 0; j < 6; j++){

				if(myL.lage == 1){  // Zopf vorne
					if(c != 3){		
						off_start = (myL.la*offLenZopf[c]/100.f);
						off_end = myL.la-(myL.la*offLenKopf[j]/100.f);
					}else{
						off_start = offLenZopf[c]; // 10 cm instead of %
						off_end = myL.la-(myL.la*offLenKopf[j]/100.f);
					}
				
					// Mittel-Linie Methode
					if(getKruemmungMittellinieLinz(&myL, off_start, off_end, diamRMI, &winkelKr)){  
						values[c][j] = myL.kr;
						degMaxKr = (int)winkelKr;
					}else{
						values[c][j] = 0;
					}

			   		if(degMaxKr < 180){
						degMaxKr = myLogs[currentLog]->degMaxKr = degMaxKr+180;
					}else{
						degMaxKr = myLogs[currentLog]->degMaxKr = 180-(360-degMaxKr);
					}


					// Mittel-Linie Methode   note: we have to invert because Zopf priori
					if(getKruemmungOberflaecheLinz(&myL, off_start, off_end, degMaxKr, diamRMI)){
						valuesMantel[c][j] = myL.kr;
					}else{
						valuesMantel[c][j] = 0;
					}
				}else if(myL.lage == 0){				//    here we have to invert  off_start && off_end for BestimmeIndex
					if(c != 3){
						off_start = myL.la-(myL.la*offLenZopf[c]/100.f);
						off_end = (myL.la*offLenKopf[j]/100.f);
					}else{
						off_start = myL.la-offLenZopf[c]; // 10 cm instead of %
						off_end = (myL.la*offLenKopf[j]/100.f);
					}
				
					// Mittel-Linie Methode
					if(getKruemmungMittellinieLinz(&myL, off_end, off_start, diamRMI, &winkelKr)){
						values[c][j] = myL.kr;
						degMaxKr = (int)winkelKr;
					}else{
						values[c][j] = 0;
					}
					

			   		if(degMaxKr < 180){
						degMaxKr = myLogs[currentLog]->degMaxKr = degMaxKr+180;
					}else{
						degMaxKr = myLogs[currentLog]->degMaxKr = 180-(360-degMaxKr);
					}



					// Mantel-Linie Methode
					if(getKruemmungOberflaecheLinz(&myL, off_end, off_start, degMaxKr, diamRMI)){
						valuesMantel[c][j] = myL.kr;
					}else{
						valuesMantel[c][j] = 0;
					}
				}else{

				}
			}
		}

		if(currentLog == 0){
			sprintf(dataSetKr.MA, "%s", S1);
		}else if(currentLog == 1){
			sprintf(dataSetKr.MA, "%s", S2);
		}

		sprintf(dataSetKr.BA, "%s", myLogs[currentLog]->mirData.ho);
		dataSetKr.tipForward = myL.lage;
		//dataSetKr.WA = 1; // this is set above
		dataSetKr.L = myLogs[currentLog]->mirData.la;
		dataSetKr.MDM = myLogs[currentLog]->mirData.rmi;
		dataSetKr.ZD = myLogs[currentLog]->mirData.zo;
		dataSetKr.A = myLogs[currentLog]->mirData.abh;
	//	dataSetKr.KR_EINSEITIG = 1;  // 1st derivative 
		dataSetKr.wKr = (int)winkelKr;

	    dataSetKr.oval = getOvalitaet(winkelKr,myLogs[currentLog]->nSlices/2);

		// Mittel-Line
		dataSetKr.mittel.PfM_3_6 = values[0][0];
		dataSetKr.mittel.PfM_3_8 = values[0][1];
		dataSetKr.mittel.PfM_3_10 = values[0][2];
		dataSetKr.mittel.PfM_3_12 = values[0][3];
		dataSetKr.mittel.PfM_3_16 = values[0][4];
		dataSetKr.mittel.PfM_3_20 = values[0][5];

		dataSetKr.mittel.PfM_5_6 = values[1][0];
		dataSetKr.mittel.PfM_5_8 = values[1][1];
		dataSetKr.mittel.PfM_5_10 = values[1][2];
		dataSetKr.mittel.PfM_5_12 = values[1][3];
		dataSetKr.mittel.PfM_5_16 = values[1][4];
		dataSetKr.mittel.PfM_5_20 = values[1][5];

		dataSetKr.mittel.PfM_8_6 = values[2][0];
		dataSetKr.mittel.PfM_8_8 = values[2][1];
		dataSetKr.mittel.PfM_8_10 = values[2][2];
		dataSetKr.mittel.PfM_8_12 = values[2][3];
		dataSetKr.mittel.PfM_8_16 = values[2][4];
		dataSetKr.mittel.PfM_8_20 = values[2][5];

		dataSetKr.mittel.PfM_10_6 = values[3][0];
		dataSetKr.mittel.PfM_10_8 = values[3][1];
		dataSetKr.mittel.PfM_10_10 = values[3][2];
		dataSetKr.mittel.PfM_10_12 = values[3][3];
		dataSetKr.mittel.PfM_10_16 = values[3][4];
		dataSetKr.mittel.PfM_10_20 = values[3][5];


		// Mantel
		dataSetKr.mantel.PfO_3_6 = valuesMantel[0][0];
		dataSetKr.mantel.PfO_3_8 = valuesMantel[0][1];
		dataSetKr.mantel.PfO_3_10 = valuesMantel[0][2];
		dataSetKr.mantel.PfO_3_12 = valuesMantel[0][3];
		dataSetKr.mantel.PfO_3_16 = valuesMantel[0][4];
		dataSetKr.mantel.PfO_3_20 = valuesMantel[0][5];

		dataSetKr.mantel.PfO_5_6 = valuesMantel[1][0];
		dataSetKr.mantel.PfO_5_8 = valuesMantel[1][1];
		dataSetKr.mantel.PfO_5_10 = valuesMantel[1][2];
		dataSetKr.mantel.PfO_5_12 = valuesMantel[1][3];
		dataSetKr.mantel.PfO_5_16 = valuesMantel[1][4];
		dataSetKr.mantel.PfO_5_20 = valuesMantel[1][5];

		dataSetKr.mantel.PfO_8_6 = valuesMantel[2][0];
		dataSetKr.mantel.PfO_8_8 = valuesMantel[2][1];
		dataSetKr.mantel.PfO_8_10 = valuesMantel[2][2];
		dataSetKr.mantel.PfO_8_12 = valuesMantel[2][3];
		dataSetKr.mantel.PfO_8_16 = valuesMantel[2][4];
		dataSetKr.mantel.PfO_8_20 = valuesMantel[2][5];

		dataSetKr.mantel.PfO_10_6 = valuesMantel[3][0];
		dataSetKr.mantel.PfO_10_8 = valuesMantel[3][1];
		dataSetKr.mantel.PfO_10_10 = valuesMantel[3][2];
		dataSetKr.mantel.PfO_10_12 = valuesMantel[3][3];
		dataSetKr.mantel.PfO_10_16 = valuesMantel[3][4];
		dataSetKr.mantel.PfO_10_20 = valuesMantel[3][5];

		logInfoDlg.setLogInfoKr(dataSetKr);
		logInfoDlg.printLogDataKr(ANALYSIS_DATA_FILE_KR);

		if(dataSetKr.WA && EXPORT_WA_PHOTOS && (currentLog == 0)){  // we export the pics of Logs having a Wurzel-Anlauf && coming form 3D1
			sprintf(img_path, "%s\\%s.jpg", LOG_PHOTO_DIR, dataSetKr.mif_file); 
			myLogs[currentLog]->pPhoto->saveImg(img_path);
		}
	}
	
}


void showLogData(bool writeToFile){
	const XRAY_SLICEINFO *xrayS;
	CString tmpStr;
	int len_;
	int div[3];

	if((myLogs[currentLog] && myLogs[logMatch]) != NULL){
		memset(&dataSet, 0, sizeof(HFAData));
		len_ = (getLengthWBark()*10/2); //myLogs[currentLog]->nSlices/2;
		dataSet.m50_len = len_;
		dataSet.m25_len = len_/2;
		dataSet.m75_len = len_+dataSet.m25_len;

		div[0] = dataSet.m50_len;
		div[1] = dataSet.m25_len;
		div[2] = dataSet.m75_len;

		int d1 = getLogDiam_WBark(areaEnd - ((areaStart-areaEnd)/2), DEGREE_MEASURE_D1);
		int d2 =  getLogDiam_WOBark(areaEnd - ((areaStart-areaEnd)/2), DEGREE_MEASURE_D1);
		int d1_2 = getLogDiam_WBark(areaEnd - ((areaStart-areaEnd)/2), DEGREE_MEASURE_D2);
		int d2_2 = getLogDiam_WOBark(areaEnd - ((areaStart-areaEnd)/2), DEGREE_MEASURE_D2);
		xrayS = getXRAY_SLICEINFO(areaEnd - ((areaStart-areaEnd)/2));

		int d_tomo, d_tomoD1, d_tomoD2;
		float d_tomoD1_woBark, d_tomoD2_woBark = 0;

		
		if(xrayS){
			dataSet.tomo.EncoderPos = xrayS->EncoderPos;
			dataSet.tomo.EncoderTicks = xrayS->EncoderTicks;
			if(myXRayLog[currentLog] != NULL){
			//	dataSet.tomo.nr_Slices = (myXRayLog[currentLog].pTomo)->getNr3DSlices();
			}
		}
		
	

		Model3D *pM = getXRayLog(); // first we perform the bark analysis
		if(pM->pTomo){
		//	(pM->pTomo)->BarkAnalysis();
		}else{
			
		}


		char *cH;
		char *cH2;
		char *date;
		char *time_;

		tmpStr = myLogs[currentLog]->logInfo_.mif_file;

		cH = tmpStr.GetBuffer(30);
		cH2 = tmpStr.GetBuffer(30);
		date = strtok(cH, "_");


		sprintf(dataSet.datum, "%s", date);
		sprintf(dataSet.uhrzeit, "%s", date);

		// MIR data
		sprintf(dataSet.mir.ho, "%s", myLogs[currentLog]->mirData.ho);
		sprintf(dataSet.mir.qk, "%s", myLogs[currentLog]->mirData.qk);
		sprintf(dataSet.mir.qs, "%s", myLogs[currentLog]->mirData.qs);
		dataSet.mir.abh = myLogs[currentLog]->mirData.abh;
		dataSet.mir.kr = myLogs[currentLog]->mirData.kr;
		dataSet.mir.la = myLogs[currentLog]->mirData.la;
		dataSet.mir.rmi = myLogs[currentLog]->mirData.rmi;
		dataSet.mir.zo = myLogs[currentLog]->mirData.zo;


		// now the actual HFA evaluation
		//############ Diam
		//3D1
		dataSet.D1.Diam_25L_M1 = getLogDiam_WBark(dataSet.m25_len, DEGREE_MEASURE_D1);
		dataSet.D1.Diam_25L_M2 = getLogDiam_WBark(dataSet.m25_len, DEGREE_MEASURE_D2);
		dataSet.D1.Diam_50L_M1 = getLogDiam_WBark(dataSet.m50_len, DEGREE_MEASURE_D1);
		dataSet.D1.Diam_50L_M2 = getLogDiam_WBark(dataSet.m50_len, DEGREE_MEASURE_D2);
		dataSet.D1.Diam_75L_M1 = getLogDiam_WBark(dataSet.m75_len, DEGREE_MEASURE_D1);
		dataSet.D1.Diam_75L_M2 = getLogDiam_WBark(dataSet.m75_len, DEGREE_MEASURE_D2);


		// 3D2
		dataSet.D2.Diam_25L_M1 = getLogDiam_WOBark(dataSet.m25_len, DEGREE_MEASURE_D1);
		dataSet.D2.Diam_25L_M2 = getLogDiam_WOBark(dataSet.m25_len, DEGREE_MEASURE_D2);
		dataSet.D2.Diam_50L_M1 = getLogDiam_WOBark(dataSet.m50_len, DEGREE_MEASURE_D1);
		dataSet.D2.Diam_50L_M2 = getLogDiam_WOBark(dataSet.m50_len, DEGREE_MEASURE_D2);
		dataSet.D2.Diam_75L_M1 = getLogDiam_WOBark(dataSet.m75_len, DEGREE_MEASURE_D1);
		dataSet.D2.Diam_75L_M2 = getLogDiam_WOBark(dataSet.m75_len, DEGREE_MEASURE_D2);

		// TOMOLOG
		dataSet.tomo.Diam_25L_M1 = getLogDiam_WBarkTOMOLOG(dataSet.m25_len, DEGREE_MEASURE_D1);
		dataSet.tomo.Diam_25L_M2 = getLogDiam_WBarkTOMOLOG(dataSet.m25_len, DEGREE_MEASURE_D2);
		dataSet.tomo.Diam_50L_M1 = getLogDiam_WBarkTOMOLOG(dataSet.m50_len, DEGREE_MEASURE_D1);
		dataSet.tomo.Diam_50L_M2 = getLogDiam_WBarkTOMOLOG(dataSet.m50_len, DEGREE_MEASURE_D2);
		dataSet.tomo.Diam_75L_M1 = getLogDiam_WBarkTOMOLOG(dataSet.m75_len, DEGREE_MEASURE_D1);
		dataSet.tomo.Diam_75L_M2 = getLogDiam_WBarkTOMOLOG(dataSet.m75_len, DEGREE_MEASURE_D2);

		//############ StdDiam +/- 5 cm
		//3D1
		int range = 50;
		dataSet.D1.StdDiam5_25L_M1 = getStdLogDiam_WBark(dataSet.m25_len,  range, DEGREE_MEASURE_D1,false);
		dataSet.D1.StdDiam5_25L_M2 = getStdLogDiam_WBark(dataSet.m25_len,  range, DEGREE_MEASURE_D2, false);
		dataSet.D1.StdDiam5_50L_M1 = getStdLogDiam_WBark(dataSet.m50_len,  range, DEGREE_MEASURE_D1, false);
		dataSet.D1.StdDiam5_50L_M2 = getStdLogDiam_WBark(dataSet.m50_len,  range, DEGREE_MEASURE_D2, false);
		dataSet.D1.StdDiam5_75L_M1 = getStdLogDiam_WBark(dataSet.m75_len,  range, DEGREE_MEASURE_D1, false);
		dataSet.D1.StdDiam5_75L_M2 = getStdLogDiam_WBark(dataSet.m75_len,  range, DEGREE_MEASURE_D2, false);

		// we take the StdDiam(5) to calc the peintinger
		dataSet.peintingerDiam = getPeintingerMeasure((int)(dataSet.D1.StdDiam5_50L_M1+dataSet.D1.StdDiam5_50L_M2)/2, dataSet.mir.ho[2]);

		// 3D2
		dataSet.D2.StdDiam5_25L_M1 = getStdLogDiam_WOBark(dataSet.m25_len,  range, DEGREE_MEASURE_D1, false);
		dataSet.D2.StdDiam5_25L_M2 = getStdLogDiam_WOBark(dataSet.m25_len,  range, DEGREE_MEASURE_D2, false);
		dataSet.D2.StdDiam5_50L_M1 = getStdLogDiam_WOBark(dataSet.m50_len,  range, DEGREE_MEASURE_D1, false);
		dataSet.D2.StdDiam5_50L_M2 = getStdLogDiam_WOBark(dataSet.m50_len,  range, DEGREE_MEASURE_D2, false);
		dataSet.D2.StdDiam5_75L_M1 = getStdLogDiam_WOBark(dataSet.m75_len,  range, DEGREE_MEASURE_D1, false);
		dataSet.D2.StdDiam5_75L_M2 = getStdLogDiam_WOBark(dataSet.m75_len,  range, DEGREE_MEASURE_D2, false);

		// TOMOLOG
		dataSet.tomo.StdDiam5_25L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m25_len, range, DEGREE_MEASURE_D1, false);
		dataSet.tomo.StdDiam5_25L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m25_len, range, DEGREE_MEASURE_D2, false);
		dataSet.tomo.StdDiam5_50L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m50_len, range, DEGREE_MEASURE_D1, false);
		dataSet.tomo.StdDiam5_50L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m50_len, range, DEGREE_MEASURE_D2, false);
		dataSet.tomo.StdDiam5_75L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m75_len, range, DEGREE_MEASURE_D1, false);
		dataSet.tomo.StdDiam5_75L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m75_len, range, DEGREE_MEASURE_D2, false);


		//############ StdDiam +/- 10 cm   here we eliminate 25% min and 25% max
		//3D1   
		range = 100;
		dataSet.D1.StdDiam10_25L_M1 = getStdLogDiam_WBark(dataSet.m25_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D1.StdDiam10_25L_M2 = getStdLogDiam_WBark(dataSet.m25_len, range, DEGREE_MEASURE_D2, true);
		dataSet.D1.StdDiam10_50L_M1 = getStdLogDiam_WBark(dataSet.m50_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D1.StdDiam10_50L_M2 = getStdLogDiam_WBark(dataSet.m50_len, range, DEGREE_MEASURE_D2, true);
		dataSet.D1.StdDiam10_75L_M1 = getStdLogDiam_WBark(dataSet.m75_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D1.StdDiam10_75L_M2 = getStdLogDiam_WBark(dataSet.m75_len, range, DEGREE_MEASURE_D2, true);

		// 3D2
		dataSet.D2.StdDiam10_25L_M1 = getStdLogDiam_WOBark(dataSet.m25_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D2.StdDiam10_25L_M2 = getStdLogDiam_WOBark(dataSet.m25_len, range, DEGREE_MEASURE_D2, true);
		dataSet.D2.StdDiam10_50L_M1 = getStdLogDiam_WOBark(dataSet.m50_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D2.StdDiam10_50L_M2 = getStdLogDiam_WOBark(dataSet.m50_len, range, DEGREE_MEASURE_D2, true);
		dataSet.D2.StdDiam10_75L_M1 = getStdLogDiam_WOBark(dataSet.m75_len, range, DEGREE_MEASURE_D1, true);
		dataSet.D2.StdDiam10_75L_M2 = getStdLogDiam_WOBark(dataSet.m75_len, range, DEGREE_MEASURE_D2,  true);

		// TOMOLOG
		dataSet.tomo.StdDiam10_25L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m25_len, range,DEGREE_MEASURE_D1, true);
		dataSet.tomo.StdDiam10_25L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m25_len, range,DEGREE_MEASURE_D2, true);
		dataSet.tomo.StdDiam10_50L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m50_len, range,DEGREE_MEASURE_D1, true);
		dataSet.tomo.StdDiam10_50L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m50_len, range,DEGREE_MEASURE_D2, true);
		dataSet.tomo.StdDiam10_75L_M1 = getStdLogDiam_WBarkTOMOLOG(dataSet.m75_len, range,DEGREE_MEASURE_D1, true);
		dataSet.tomo.StdDiam10_75L_M2 = getStdLogDiam_WBarkTOMOLOG(dataSet.m75_len, range,DEGREE_MEASURE_D2, true);

		
		

			// TOMOLOG Bark Detection  we take the data form 3 differen positions
			//### L*0.50
				// MuMax

				dataSet.tomo.pos50.D_MuMax_1 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5);
				dataSet.tomo.pos50.D_MuMax_2 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5);

				// MuD
				dataSet.tomo.pos50.D_MuD_1 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1);
				dataSet.tomo.pos50.D_MuD_2 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1);

				// MuD2
				dataSet.tomo.pos50.D_MuD2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2);
				dataSet.tomo.pos50.D_MuD2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2);

				// ST
				dataSet.tomo.pos50.D_ST_1 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3);
				dataSet.tomo.pos50.D_ST_2 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3);

				// ST2
				dataSet.tomo.pos50.D_ST2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4);
				dataSet.tomo.pos50.D_ST2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[0], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4);


			   //### L*0.25
				// MuMax
				dataSet.tomo.pos25.D_MuMax_1 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5);
				dataSet.tomo.pos25.D_MuMax_2 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5);

				// MuD
				dataSet.tomo.pos25.D_MuD_1 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1);
				dataSet.tomo.pos25.D_MuD_2 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1);

				// MuD2
				dataSet.tomo.pos25.D_MuD2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2);
				dataSet.tomo.pos25.D_MuD2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2);

				// ST
				dataSet.tomo.pos25.D_ST_1 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3);
				dataSet.tomo.pos25.D_ST_2 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3);

				// ST2
				dataSet.tomo.pos25.D_ST2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4);
				dataSet.tomo.pos25.D_ST2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[1], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4);


				//### L*0.75
				dataSet.tomo.pos75.D_MuMax_1 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5);
				dataSet.tomo.pos75.D_MuMax_2 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5);

				// MuD
				dataSet.tomo.pos75.D_MuD_1 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1);
				dataSet.tomo.pos75.D_MuD_2 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1);

				// MuD2
				dataSet.tomo.pos75.D_MuD2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2);
				dataSet.tomo.pos75.D_MuD2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2);

				// ST
				dataSet.tomo.pos75.D_ST_1 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3);
				dataSet.tomo.pos75.D_ST_2 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3);

				// ST2
				dataSet.tomo.pos75.D_ST2_1 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4);
				dataSet.tomo.pos75.D_ST2_2 = (int)getLogDiam_WOBarkTOMOLOG(div[2], DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4);


					
			
			// TOMOLOG   range +/- 10 cm
			// L*0.50   M1
			range = 50;
			dataSet.tomo.pos50.MW5_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos50.MW5_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos50.MW5_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos50.MW5_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos50.MW5_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, false);  // ST2


			// L*0.50   M12
			dataSet.tomo.pos50.MW5_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos50.MW5_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos50.MW5_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos50.MW5_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos50.MW5_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, false);  // ST2


			// L*0.25   M1

			dataSet.tomo.pos25.MW5_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos25.MW5_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos25.MW5_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos25.MW5_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos25.MW5_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, false);  // ST2


			// L*0.25   M12
			dataSet.tomo.pos25.MW5_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos25.MW5_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos25.MW5_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos25.MW5_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos25.MW5_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, false);  // ST2
			

			// L*0.75   M1

			dataSet.tomo.pos75.MW5_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos75.MW5_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos75.MW5_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos75.MW5_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos75.MW5_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, false);  // ST2


			// L*0.75   M12
			dataSet.tomo.pos75.MW5_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, false);  //MuMax
			dataSet.tomo.pos75.MW5_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, false);  // MuD
			dataSet.tomo.pos75.MW5_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, false);  // MuD2
			dataSet.tomo.pos75.MW5_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, false);  // ST
			dataSet.tomo.pos75.MW5_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, false);  // ST2
			


			//#################################################################
			// TOMOLOG range +/- 10 cm
			range = 100;
			dataSet.tomo.pos50.MW10_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos50.MW10_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos50.MW10_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos50.MW10_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos50.MW10_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, true);  // ST2


			// L*0.50   M2
			dataSet.tomo.pos50.MW10_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos50.MW10_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos50.MW10_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos50.MW10_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos50.MW10_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[0], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, true);  // ST2


			// L*0.25   M1

			dataSet.tomo.pos25.MW10_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos25.MW10_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos25.MW10_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos25.MW10_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos25.MW10_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, true);  // ST2


			// L*0.25   M2
			dataSet.tomo.pos25.MW10_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos25.MW10_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos25.MW10_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos25.MW10_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos25.MW10_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[1], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, true);  // ST2
			

			// L*0.75   M1

			dataSet.tomo.pos75.MW10_MuMax_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos75.MW10_MuD_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos75.MW10_MuD2_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos75.MW10_ST_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos75.MW10_ST2_1 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4, true);  // ST2


			// L*0.75   M2
			dataSet.tomo.pos75.MW10_MuMax_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5, true);  //MuMax
			dataSet.tomo.pos75.MW10_MuD_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1, true);  // MuD
			dataSet.tomo.pos75.MW10_MuD2_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2, true);  // MuD2
			dataSet.tomo.pos75.MW10_ST_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3, true);  // ST
			dataSet.tomo.pos75.MW10_ST2_2 = getStdLogDiam_WOBarkTOMOLOG(div[2], range, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4, true);  // ST2
			
			populateONORM(dataSet);
		
		


		
		//#############################################


		myLogs[currentLog]->logInfo_.D1_diam1 = d1;
		myLogs[currentLog]->logInfo_.D2_diam1 = d2;
		myLogs[currentLog]->logInfo_.D1_diam2 = d1_2;
		myLogs[currentLog]->logInfo_.D2_diam2 = d2_2;
		myLogs[currentLog]->logInfo_.len_wBark = getLengthWBark();
		myLogs[currentLog]->logInfo_.len_woBark = getLengthWOBark();
		myLogs[currentLog]->logInfo_.TOMOLOG_diam = d_tomo;
		myLogs[currentLog]->logInfo_.TOMOLOG_diamD1 = d_tomoD1;
		myLogs[currentLog]->logInfo_.TOMOLOG_diamD2 = d_tomoD2;
		myLogs[currentLog]->logInfo_.TOMOLOG_D1_woBark = d_tomoD1_woBark;
		myLogs[currentLog]->logInfo_.TOMOLOG_D2_woBark = d_tomoD2_woBark;

		logInfoDlg.setLogInfo(myLogs[currentLog]->logInfo_, dataSet);
		
		if(!writeToFile){
			logInfoDlg.DoModal(); 
		}else{
			logInfoDlg.printLogData(ANALYSIS_DATA_FILE);
		}
	}
 }


 void populateONORM(HFAData &data){
	int md_s1;
	int md_s2;
	int md_tomo;
	int md_MuMax;
	int md_MuD;
	int md_MuD2;
	int md_ST;
	int md_ST2;

	int off_md_s1;
	int off_md_s2;
	int off_md_tomo;
	int off_md_MuMax;
	int off_md_MuD;
	int off_md_MuD2;
	int off_md_ST;
	int off_md_ST2; 

	int start_off = data.m50_len-150;  // +/- 15 cm
	int end_off = data.m50_len+150;

	md_s1 = getMinimumDiamDiShapeWBark(start_off, end_off, off_md_s1);
	md_s2 = getMinimumDiamDiShapeWOBark(start_off, end_off, off_md_s2);
	md_tomo = getMinimumDiamTOMO(start_off, end_off, off_md_tomo, (XRAY_DIAMTYPE)0);
	md_MuMax =  getMinimumDiamTOMO(start_off, end_off, off_md_MuMax, (XRAY_DIAMTYPE)5);
	md_MuD = getMinimumDiamTOMO(start_off, end_off, off_md_MuD, (XRAY_DIAMTYPE)1);
	md_MuD2 = getMinimumDiamTOMO(start_off, end_off, off_md_MuD2, (XRAY_DIAMTYPE)2); 
	md_ST = getMinimumDiamTOMO(start_off, end_off, off_md_ST, (XRAY_DIAMTYPE)3); 
	md_ST2 = getMinimumDiamTOMO(start_off, end_off, off_md_ST2, (XRAY_DIAMTYPE)4); 

	// we store the diameters
	data.norm.MD_S1 = md_s1;
	data.norm.MD_S2 = md_s2;
	data.norm.MD_TOMO = md_tomo;
	data.norm.MD_MuMax = md_MuMax;
	data.norm.MD_MuD = md_MuD;
	data.norm.MD_MuD2 = md_MuD2;
	data.norm.MD_ST = md_ST;
	data.norm.MD_ST2 = md_ST2;

	// the offsets where we found the min
	data.norm.L_MD_S1 = off_md_s1;
	data.norm.L_MD_S2 = off_md_s2;
	data.norm.L_MD_TOMO = off_md_tomo;
	data.norm.L_MD_MuMax = off_md_MuMax;
	data.norm.L_MD_MuD = off_md_MuD;
	data.norm.L_MD_MuD2 = off_md_MuD2;
	data.norm.L_MD_ST = off_md_ST;
	data.norm.L_MD_ST2 = off_md_ST2;



	// the single alg diameters for the offsets below
	data.norm.D_LMD_S1_1 = getLogDiam_WBark(data.m50_len+off_md_s1, DEGREE_MEASURE_D1);
	data.norm.D_LMD_S2_1 = getLogDiam_WOBark(data.m50_len+off_md_s2, DEGREE_MEASURE_D1);
	data.norm.D_LMD_TOMO_1 = getLogDiam_WBarkTOMOLOG(data.m50_len+off_md_tomo, DEGREE_MEASURE_D1);		
	data.norm.D_LMD_S1_2 = getLogDiam_WBark(data.m50_len+off_md_s1, DEGREE_MEASURE_D2);
	data.norm.D_LMD_S2_2 = getLogDiam_WOBark(data.m50_len+off_md_s2, DEGREE_MEASURE_D2);
	data.norm.D_LMD_TOMO_2 = getLogDiam_WBarkTOMOLOG(data.m50_len+off_md_tomo, DEGREE_MEASURE_D2);

	// Algs
	data.norm.D_LMD_MuD_1 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuD, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)1); 
	data.norm.D_LMD_MuD2_1 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuD2, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)2); 
	data.norm.D_LMD_ST_1 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_ST, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)3); 
	data.norm.D_LMD_ST2_1 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_ST2, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)4); 
	data.norm.D_LMD_MuMax_1 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuMax, DEGREE_MEASURE_D1, (XRAY_DIAMTYPE)5); 

	data.norm.D_LMD_MuD_2 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuD, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)1); 
	data.norm.D_LMD_MuD2_2 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuD2, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)2); 
	data.norm.D_LMD_ST_2 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_ST, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)3); 
	data.norm.D_LMD_ST2_2 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_ST2, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)4); 
	data.norm.D_LMD_MuMax_2 = getLogDiam_WOBarkTOMOLOG(data.m50_len+off_md_MuMax, DEGREE_MEASURE_D2, (XRAY_DIAMTYPE)5); 


 }

 int getMinimumDiamTOMO(int start, int end, int &off_ret, XRAY_DIAMTYPE type){
	const XRAY_SLICEINFO *xrayS;
	const int size__ = (end-start);
	static int map[310][2];
	static int val[310];
	int minDiam = 0;
	int diam1, diam2;
	int off = 0;


	// we loop +/- 15 cm arount the middle
	for(int i = start; i <= end; i++){ 
		    xrayS = getXRAY_SLICEINFO(i);  // 50*l
			if(xrayS){
				switch(type){
				   case Optic:
					    diam1 =  (int)xrayS->RSx[0] + xrayS->RDx[0]*10;
						diam2 =  (int)xrayS->RSx[1] + xrayS->RDx[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;
				   case MuD:
					    diam1 =  (int)xrayS->RSxMuD[0] + xrayS->RDxMuD[0]*10;
						diam2 =  (int)xrayS->RSxMuD[1] + xrayS->RDxMuD[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;
				   case MuD2:
					    diam1 =  (int)xrayS->RSxMuD2[0] + xrayS->RDxMuD2[0]*10;
						diam2 =  (int)xrayS->RSxMuD2[1] + xrayS->RDxMuD2[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;
				   case ST:
					    diam1 =  (int)xrayS->RSxST[0] + xrayS->RDxST[0]*10;
						diam2 =  (int)xrayS->RSxST[1] + xrayS->RDxST[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;
				   case ST2:
					    diam1 =  (int)xrayS->RSxST2[0] + xrayS->RDxST2[0]*10;
						diam2 =  (int)xrayS->RSxST2[1] + xrayS->RDxST2[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;
				   case MuMax:
					    diam1 =  (int)xrayS->RSxMuMax[0] + xrayS->RDxMuMax[0]*10;
						diam2 =  (int)xrayS->RSxMuMax[1] + xrayS->RDxMuMax[1]*10;
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;
					   break;


				}

			}
			xrayS = NULL;

	}
	qsort(val, size__, sizeof(int), compareMin);
	minDiam = val[0];  // this is now the minimum
	
	// now get the offset where we found the minimum
	for(int x = 0; x <= size__; x++){
		if(map[x][0] == minDiam){
			off_ret = map[x][1];
			if(x < 150){
				off_ret *= -1;
			}else if(x >= 150){
				off_ret = off_ret%150;

			}
			break;
		}
	}
	return minDiam;
 }


  int getMinimumDiamDiShapeWBark(int start, int end, int &off_ret){
	const XRAY_SLICEINFO *xrayS;
	const int size__ = (end-start);
	static int map[310][2];
	static int val[310];
	int minDiam = 0;
	int diam1, diam2;
	int off = 0;

	// we loop +/- 15 cm arount the middle
	for(int i = start; i <= end; i++){ 
			if((myLogs[currentLog] != NULL)){
					    diam1 =  (int)myLogs[currentLog]->pData->GetDiameterAt(i, DEGREE_MEASURE_D1, bClambDiShapeDiam);
						diam2 =  (int)myLogs[currentLog]->pData->GetDiameterAt(i, DEGREE_MEASURE_D2, bClambDiShapeDiam);
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = i;
						val[off] = map[off][0];
						off++;

			}
			xrayS = NULL;

	}
	qsort(val, size__, sizeof(int), compareMin);
	minDiam = val[0];  // this is now the minimum
	
	// now get the offset where we found the minimum
	for(int x = 0; x <= size__; x++){
		if(map[x][0] == minDiam){
			off_ret = map[x][1];
			if(x < 150){
				off_ret *= -1;
			}else if(x >= 150){
				off_ret = off_ret%150;

			}
			break;
		}
	}
	return minDiam;
 }

  int getMinimumDiamDiShapeWOBark(int start, int end, int &off_ret){
	const XRAY_SLICEINFO *xrayS;
	const int size__ = (end-start);
	static int map[310][2];
	static int val[310];
	int minDiam = 0;
	int diam1, diam2;
	int off = 0;

	// we loop +/- 15 cm arount the middle
	for(int i = start; i <= end; i++){ 
			if((myLogs[logMatch] != NULL)){
					    diam1 =  (int)myLogs[logMatch]->pData->GetDiameterAt(i, DEGREE_MEASURE_D1, bClambDiShapeDiam);
						diam2 =  (int)myLogs[logMatch]->pData->GetDiameterAt(i, DEGREE_MEASURE_D2, bClambDiShapeDiam);
					   	map[off][0] = (diam1+diam2)/2;
						map[off][1] = off;
						val[off] = map[off][0];
						off++;

			}
			xrayS = NULL;

	}
	qsort(val, size__, sizeof(int), compareMin);
	minDiam = val[0];  // this is now the minimum
	
	// now get the offset where we found the minimum
	for(int x = 0; x <= size__; x++){
		if(map[x][0] == minDiam){
			off_ret = map[x][1];
			if(x < 150){
				off_ret *= -1;
			}else if(x >= 150){
				off_ret = off_ret%150;

			}
			break;
		}
	}
	return minDiam;
 }

 int compareMin(const void * a, const void * b){
	 return ( *(int*)a - *(int*)b );
 }

 int getPeintingerMeasure(int diam, char type){
   int dor, abz = 0;

   
   switch(type){
		   case 'F':
				dor = (FAKTOR_FICHTE*diam + OFFSET_FICHTE + 50000)/100000;
			   break;
		   case 'L':
				 dor = (FAKTOR_LAERCHE*diam + OFFSET_LAERCHE + 50000)/100000;
			   break;
		   case 'K':
				dor = (FAKTOR_KIEFER*diam + OFFSET_KIEFER + 50000)/100000;
			   break;
		   case 'T':
				dor = (FAKTOR_KIEFER*diam + OFFSET_KIEFER + 50000)/100000;
			   break;
   }
   
  
   return dor;
 }


 /*
		this function is used within qsort(..) for Kr calculation
  */
int ConfAbst(const void *e1,const void *e2){
   int *a,*b;

   a=(int*)e1;
   b=(int*)e2;

   if ((*a)<(*b)) return 1;
   if ((*a)>(*b)) return -1;
   return 0;
}


bool isMeasureOverrun(const XRAY_SLICEINFO *pS){
	int det1sx = 0;
	int det1dx = 0;
	int det2sx = 0;
	int det2dx = 0;

	det1sx = TOTAL_NR_PIXELS-pS->Sx[0];
	det1dx = TOTAL_NR_PIXELS-pS->Dx[0];

	det2sx = TOTAL_NR_PIXELS-pS->Sx[1];
	det2dx = TOTAL_NR_PIXELS-pS->Dx[1];

	if((det1sx > MIN_PIXEL_REST) && (det1dx > MIN_PIXEL_REST) && (det2sx > MIN_PIXEL_REST) && (det2dx > MIN_PIXEL_REST)  ){
		return false;
	}else{
		return true;
	}
}



void VisualizeKruemmung(){
	static char m__[30];
	static int offLenZopf[100];
	static int offLenKopf[100];
	static int values[25][25];
	static int valuesMantel[25][25];
	static char img_path[256];
	double winkelKr = 0;
	int einseitigeKr;
	int cntWA = 0;
	int krTmp, krTmp2, krTmp3, krTmp4;
	CString S1 = "S1";
	CString S2 = "S2";

	int off_start = 0;
	int off_end = 0;
	memset(&dataSetKr, 0, sizeof(KRData));
	memset(&offLenZopf, 0, sizeof(offLenZopf));
	memset(&offLenKopf, 0, sizeof(offLenKopf));
	memset(&myL, 0, sizeof(LOG));
	memset(&values, 0, sizeof(values));
	memset(&valuesMantel, 0, sizeof(valuesMantel));

	int degMaxKr = 0;
	// len offsets in % from Zopf
	offLenZopf[0] = 3;
	offLenZopf[1] = 5;
	offLenZopf[2] = 8;
	offLenZopf[3] = 100; // this is indeed 10 cm  -> 
				
	// len offsets in % from Kopf
	offLenKopf[0] = 6;
	offLenKopf[1] = 8;
	offLenKopf[2] = 10;
	offLenKopf[3] = 12;
	offLenKopf[4] = 16;
	offLenKopf[5] = 20;

	if(myLogs[currentLog] != NULL && myLogs[logMatch] != NULL){
		myL.la = myLogs[currentLog]->getLength()*10;  // we procceed in mm
		myL.nrslices = myLogs[currentLog]->nSlices;
		myL.lage =   (myLogs[0]->pData)->IsTipForward(); // Zopf priori   here we always take 3D1 as measure  myLogs[0]

		int diamRMI = myLogs[currentLog]->mirData.rmi; 


	//	smoothLog(myLogs[currentLog]);   // we smooth the surface
		for(int i = 0;i < myLogs[currentLog]->nSlices; i++){
			myL.slices[i].centerx = myLogs[currentLog]->Mess3D[i].zx;
			myL.slices[i].centery = myLogs[currentLog]->Mess3D[i].zy;
			myL.slices[i].z = myLogs[currentLog]->Mess3D[i].x;

			for(int j =0; j < 360; j++){
				myL.slices[i].fpoint[j].x = myLogs[currentLog]->Mess3D[i].pointx[j];
				myL.slices[i].fpoint[j].y = myLogs[currentLog]->Mess3D[i].pointy[j];

			}

			/*
			int sx = 0;
			int sy = 0;
			for(int z =0; z < 360; z++){
				myL.slices[i].fpoint[z].x = myLogs[currentLog]->Mess3D[i].pointx[z];
				myL.slices[i].fpoint[z].y = myLogs[currentLog]->Mess3D[i].pointy[z];
				sx+=myL.slices[i].fpoint[z].x;
				sy+=myL.slices[i].fpoint[z].y;
			}
			myL.slices[i].centerx=sx/360;			// we also need to adjust the central point
			myL.slices[i].centery=sy/360;
			*/
			
		}




		dataSetKr.KR_EINSEITIG = hasKrEinseitig(&myL, 10, myL.nrslices-10);

		// we reset the offsets
		off_end = 0;
		off_start = 0;
		
		myL.einseitigeKr_90  = 0;
		myL.einseitigeKr_180 = 0;
		
		// this algorithm is used to detect a "Einseitiger-Wurzelalnauf"		
		int steigungB, steigungD, begX, begY, endX, endY, ptA, ptB = 0;
		int invertedIndex_1, invertedIndex_2 = 0;
		if(myL.lage){ // Zopf at beginning
			invertedIndex_1 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF;
			invertedIndex_2 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_B;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery
				;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungB = ptA/ptB;
			}else{
				steigungB = ptA/1;
			}

			invertedIndex_1 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_C;
			invertedIndex_2 = (myLogs[currentLog]->nSlices)-OFFSET_WURZEL_ANLAUF_D;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungD = ptA/ptB;
			}else{
				steigungD = ptA/1;
			}

			if(steigungB > steigungD+DIFF_OFFSET_WURZEL_ANLAUF){
				dataSetKr.WA = 1;
			}else{
				dataSetKr.WA = 0;
			}
		}else{	// Zopf at end
			invertedIndex_1 = OFFSET_WURZEL_ANLAUF;
			invertedIndex_2 = OFFSET_WURZEL_ANLAUF_B;
			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungB = ptA/ptB;
			}else{
				steigungB = ptA/1;
			}

			invertedIndex_1 = OFFSET_WURZEL_ANLAUF_C;
			invertedIndex_2 = OFFSET_WURZEL_ANLAUF_D;

			begX = myL.slices[invertedIndex_1].centerx;
			begY = myL.slices[invertedIndex_1].centery;
			endX = myL.slices[invertedIndex_2].centerx;
			endY = myL.slices[invertedIndex_2].centery;

			ptA = abs(begX-endX);
			ptB = abs(begY-endY);

			if(ptB != 0){
				steigungD = ptA/ptB;
			}else{
				steigungD = ptA/1;
			}

			if(steigungB > steigungD+DIFF_OFFSET_WURZEL_ANLAUF){
				dataSetKr.WA = 1;
			}else{
				dataSetKr.WA = 0;
			}

		}
		// we store both WA-factors in the Dataset
		dataSetKr.wurzelAn.stB = steigungB;
		dataSetKr.wurzelAn.stD = steigungD;
		dataSetKr.wurzelAn.stDiff = steigungB-steigungD;  // if negative -> steigungD > steigungB

		for(int c = 0; c < 4; c++){
			for(int j = 0; j < 6; j++){
				if(myL.lage == 1){  // Zopf vorne
					if(c != 3){		
						off_start = (myL.la*offLenZopf[c]/100.f);
						off_end = myL.la-(myL.la*offLenKopf[j]/100.f);
					}else{
						off_start = offLenZopf[c]; // 10 cm instead of %
						off_end = myL.la-(myL.la*offLenKopf[j]/100.f);
					}
				
					// Mittel-Linie Methode
					if(getKruemmungMittellinieLinz(&myL, off_start, off_end, diamRMI, &winkelKr)){  
						values[c][j] = myL.kr;
						degMaxKr = (int)winkelKr;
		
					}else{
						values[c][j] = 0;
					}

					
			   		if(degMaxKr < 180){
						degMaxKr = degMaxKr+180;
					}else{
						degMaxKr =  180-(360-degMaxKr);
					}




					// Mittel-Linie Methode   note: we have to invert because Zopf priori
					if(getKruemmungOberflaecheLinz(&myL, off_start, off_end, degMaxKr, diamRMI)){
						valuesMantel[c][j] = myL.kr;
					}else{
						valuesMantel[c][j] = 0;
					}
				}else if(myL.lage == 0){				//    here we have to invert  off_start && off_end for BestimmeIndex
					if(c != 3){
						off_start = myL.la-(myL.la*offLenZopf[c]/100.f);
						off_end = (myL.la*offLenKopf[j]/100.f);
					}else{
						off_start = myL.la-offLenZopf[c]; // 10 cm instead of %
						off_end = (myL.la*offLenKopf[j]/100.f);
					}
				
					// Mittel-Linie Methode
					if(getKruemmungMittellinieLinz(&myL, off_end, off_start, diamRMI, &winkelKr)){
						values[c][j] = myL.kr;
						degMaxKr = (int)winkelKr;
		
					}else{
						values[c][j] = 0;
					}	

					if(degMaxKr < 180){
						degMaxKr = degMaxKr+180;
					}else{
						degMaxKr =  180-(360-degMaxKr);
					}


					// Mantel-Linie Methode
					if(getKruemmungOberflaecheLinz(&myL, off_end, off_start, degMaxKr, diamRMI)){
						valuesMantel[c][j] = myL.kr;
					}else{
						valuesMantel[c][j] = 0;
					}
				}else{

				}

			}
		}


			   myLogs[logMatch]->degMaxKr = degMaxKr;
			   sprintf(m__, "%d", degMaxKr); 
			   MSG(m__);

/*
			   	if(degMaxKr < 180){
							myLogs[currentLog]->degMaxKr = degMaxKr+180;
				}else{
					myLogs[currentLog]->degMaxKr = 180-(360-degMaxKr);
				}
				*/

				myLogs[currentLog]->degMaxKr = degMaxKr;
			    myLogs[currentLog]->lineKr.krAvail = true;
				
				renderLog(1);  
				renderLog(2);

	}		


}

 BOOL getKruemmungMittellinieBx(LOG *log, int off_start, int off_end){
					int x,y,xr,yr,x1,y1,z1,x2,y2,z2,z;
					int iv,ib,lv,lb;
					int i,nr,d;
					int abst[MAX_SLICE];

					log->kr=0;


						lv=(off_start);
						lb=(off_start+100);
	
					BestimmeIndex(log,lv,lb,&iv,&ib);

					nr=0;
					x1=0;
					y1=0;
					for (i=iv; i<=ib; i++) {
						x1+=log->slices[i].centerx;
						y1+=log->slices[i].centery;
						nr++;
					}
					if (nr>0) {
						x1/=nr;
						y1/=nr;
					} else {
						return FALSE;
					}

					z1=(lv+lb)/2;

						lv=off_end-100;
						lb=off_end;

						
						BestimmeIndex(log,lv,lb,&iv,&ib);

					nr=0;
					x2=0;
					y2=0;
					for (i=iv; i<=ib; i++) {
						x2+=log->slices[i].centerx;
						y2+=log->slices[i].centery;
						nr++;
					}
					if (nr>0) {
						x2/=nr;
						y2/=nr;
					} else {
						return FALSE;
					}

					z2=(lv+lb)/2;
					
					nr=0;
					for (i=0; i<log->nrslices; i++) {
						if (log->slices[i].z>z1 && log->slices[i].z<z2) {
							xr=log->slices[i].centerx;
							yr=log->slices[i].centery;
							z=log->slices[i].z;
							x=((long)((z-z1)*(x2-x1)*100l)/(z2-z1)+x1*100l)/100l;
							y=((long)((z-z1)*(y2-y1)*100l)/(z2-z1)+y1*100l)/100l;   
							d=(x-xr)*(x-xr)+(y-yr)*(y-yr);  
							abst[nr]=d;
							nr++;
						}

					}
   
					qsort(abst,nr,sizeof(int),ConfAbst);
   
					d=0;
					nr=0;
					for (i=3; i<10; i++) {
					   d+=abst[i];
					   nr++;
				   }

					if (nr>0){ d/=nr;   }
                                 
				   log->kr=(int)((sqrt((double)d)/10.0+0.5)*10);  
				   log->kr*=100l;
				   if (log->la!=0) log->kr/=log->la;
				   else log->kr=0;

				   

				   return TRUE;
 }

  BOOL getKruemmungMittellinieLinz(LOG *log, int off_start, int off_end, int diam__, double *dWinkel){
 	double x,y,z,xr,yr,x1,y1,z1,x2,y2,z2,oldz,dx,dy, zz1, zz2;
	int iv,ib,lv,lb;

	int i,k,nr;

	double eta, etaMin, etaCurv, etaFinal, etaTurn;
	double s, s1, mins;
	double Ca, Sa;
	double kruemmung;
	double d2, maxd2 = 0;


	// Mittelpunkt Anfang


	int avgPt = 0;
	avgPt = POS_AVG_OFFSET/2;

	lv= off_start-avgPt;
	lb = off_start+avgPt;


	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x1=0;
	y1=0;
	for (i=iv; i<=ib; i++) {
		x1+=(double)log->slices[i].centerx;
		y1+=(double)log->slices[i].centery;

		nr++;
	}
	zz1 = (double)log->slices[iv].z;

	if (nr>0) {
		x1/=(double)nr;
		y1/=(double)nr;
	}
	else {
		return FALSE;
	}
	z1=(double)(lv+lb)/2.;


	lv= off_end-avgPt;
	lb = off_end+avgPt;

	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x2=0;
	y2=0;
	for (i=iv; i<=ib; i++) {
		x2+=(double)log->slices[i].centerx;
		y2+=(double)log->slices[i].centery;
		nr++;
	}

	zz2 = (double)log->slices[iv].z;

	if (nr>0) {
		x2/=(double)nr;
		y2/=(double)nr;
	}
	else {
		return FALSE;
	}
	z2=(double)(lv+lb)/2;


	log->rx[0]=(int)x1;
	log->ry[0]=(int)y1;
	log->rz[0]=(int)z1;

	log->rx[1]=(int)x2;
	log->ry[1]=(int)y2;
	log->rz[1]=(int)z2;

	if(!myLogs[currentLog]->lineKr.krAvail){
		// we store the Segments for visualization purposes
		myLogs[currentLog]->lineKr.m_1_x = ((float)x1)/10000.0f;
		myLogs[currentLog]->lineKr.m_1_y = ((float)y1)/10000.0f;
		myLogs[currentLog]->lineKr.m_1_z = ((float)zz1)/10000.0f;
		myLogs[currentLog]->lineKr.m_2_x = ((float)x2)/10000.0f;
		myLogs[currentLog]->lineKr.m_2_y = ((float)y2)/10000.0f;
		myLogs[currentLog]->lineKr.m_2_z = ((float)zz2)/10000.0f;
	}

	etaMin=0;

	mins=1e10; 

	for (eta=0.001; eta<=1.58; eta+=0.025) {
		for (k=0; k<2; k++) {
			Sa=sin(eta);    
			if (k==1) Sa=-Sa;
			Ca=cos(eta);
			s=0;  
			oldz=z1;   

			nr=0;
			for (i=0; i<log->nrslices; i++) {
				if (log->slices[i].z>z1 && log->slices[i].z<z2) {

					xr=log->slices[i].centerx;
					yr=log->slices[i].centery;
					z=log->slices[i].z;

					x=x1+((z-z1)*(x2-x1))/(z2-z1);
					y=y1+((z-z1)*(y2-y1))/(z2-z1);

					dx=(xr-x)*Ca-(yr-y)*Sa;
					dy=(xr-x)*Sa+(yr-y)*Ca;

					d2 = dx*dx+dy*dy;

					if(d2 > maxd2){
						maxd2 = d2;
					}

					s+=((double)abs((int)dx)*(z-oldz));

					nr++;

					oldz=z;
				}
			}

			if (oldz>z1) {
				s/=(oldz-z1);
			}
			else {
				s=0;
			}



			if (s>0 && s<mins) {
				mins=s;
				if (k==0) etaMin=eta;
				else etaMin=-eta;
			}
		}
	}

	
	
	if(maxd2 > 0 && (z2 -z1 > 0)){
		kruemmung = sqrt(maxd2); 
	}else{
		kruemmung = 0;
	}


		// now we calc the Angle
	eta=etaMin;

	Sa=sin(eta);    
	Ca=cos(eta);
	s=0;  
	s1=0;

	oldz=z1;

	for (i=0; i<log->nrslices; i++) {
		if (log->slices[i].z>z1 && log->slices[i].z<z2) {

			xr=log->slices[i].centerx;
			yr=log->slices[i].centery;
			z=log->slices[i].z;

			x=x1+((z-z1)*(x2-x1))/(z2-z1);
			y=y1+((z-z1)*(y2-y1))/(z2-z1);

			dx=(xr-x)*Ca-(yr-y)*Sa;
			dy=(xr-x)*Sa+(yr-y)*Ca;

			if (dy>0) {
				s+=dy*(z-oldz);   
			}
			else if (dy<0) {                
				s1+=(-dy)*(z-oldz);    
			}  
			oldz=z;
		}

	}       


	if (s1>s) {  
		if (etaMin>0) {
			etaMin-=3.1415;
		}
		else {
			etaMin+=3.1415;
		}
	}

	etaCurv=3.1415/2.0-etaMin;				// ermittelter Kruemmungswinkel
	if (etaCurv<-3.1415) etaCurv+=2.0*3.1415;
	else if (etaCurv>3.1415) etaCurv-=2.0*3.1415;

	etaFinal=(double)90*3.1415/180.;		// gewuenschter Kruemmungswinkel
	if (etaFinal<-3.1415) etaFinal+=2.0*3.1415;
	else if (etaFinal>3.1415) etaFinal-=2.0*3.1415;

	*dWinkel = abs((int)(etaCurv*(180/PI))); //Uhr2Grad(radMath2GradUhr(etaCurv));




    diam__ *= 10;		// 1/10 mm

	if(USE_PERCENTAGE_DIAM_RATION){
		log->kr = (int)((100*(kruemmung))/diam__);  // in % ratio of the standard diameter
	}else{
		log->kr = (int)kruemmung;  // in 1/10 mm
	}
	
	return TRUE;
  }


  /*
		Krümmung nach der Mantel-Linie
  */
  BOOL getKruemmungOberflaecheLinz(LOG *log, int off_start, int off_end, int degMaxKr, int diam__){
 	double x,y,z,xr,yr,x1,y1,z1,x2,y2,z2,oldz,dx,dy;
	int iv,ib,lv,lb;

	int i,k,nr;

	double eta, etaMin, etaCurv, etaFinal, etaTurn;
	double s, s1, mins;
	double Ca, Sa;
	double kruemmung;
	double d2, maxd2 = 0;


	int avgPt = 0;
	avgPt = POS_AVG_OFFSET/2;

	lv= off_start-avgPt;
	lb = off_start+avgPt;


	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x1=0;
	y1=0;
	for (i=iv; i<=ib; i++) {
		x1+=log->slices[i].fpoint[degMaxKr].x;
		y1+=log->slices[i].fpoint[degMaxKr].y;
		nr++;
	}
	if (nr>0) {
		x1/=(double)nr;
		y1/=(double)nr;
	}
	else {
		return FALSE;
	}
	z1=(double)(lv+lb)/2.;


	// Mittelpunkt Ende

	lv= off_end-avgPt;
	lb = off_end+avgPt;

	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x2=0;
	y2=0;
	for (i=iv; i<=ib; i++) {
		x2+=log->slices[i].fpoint[degMaxKr].x;
		y2+=log->slices[i].fpoint[degMaxKr].y;
		nr++;
	}
	if (nr>0) {
		x2/=(double)nr;
		y2/=(double)nr;
	}
	else {
		return FALSE;
	}
	z2=(double)(lv+lb)/2;


	log->rx[0]=(int)x1;
	log->ry[0]=(int)y1;
	log->rz[0]=(int)z1;

	log->rx[1]=(int)x2;
	log->ry[1]=(int)y2;
	log->rz[1]=(int)z2;


	if(!myLogs[currentLog]->lineKr.krAvail){
		// we store the Segments for visualization purposes
		myLogs[currentLog]->lineKr.o_1_x = (float)x1/10000.0f;
		myLogs[currentLog]->lineKr.o_1_y = (float)y1/10000.0f;
		myLogs[currentLog]->lineKr.o_1_z = (float)z1/10000.0f;
		myLogs[currentLog]->lineKr.o_2_x = (float)x2/10000.0f;
		myLogs[currentLog]->lineKr.o_2_y = (float)y2/10000.0f;
		myLogs[currentLog]->lineKr.o_2_z = (float)z2/10000.0f;
		myLogs[currentLog]->lineKr.krAvail = true;
	}



	etaMin=0;

	mins=1e10; 

	
	for (eta=0.001; eta<=1.58; eta+=0.025) {
		for (k=0; k<2; k++) {
			Sa=sin(eta);    
			if (k==1) Sa=-Sa;
			Ca=cos(eta);
			


			s=0;  
			oldz=z1;   

			nr=0;
			for (i=0; i<log->nrslices; i++) {
				if (log->slices[i].z>z1 && log->slices[i].z<z2) {

					xr=log->slices[i].fpoint[degMaxKr].x;
					yr=log->slices[i].fpoint[degMaxKr].y;
					z=log->slices[i].z;

					x=x1+((z-z1)*(x2-x1))/(z2-z1);
					y=y1+((z-z1)*(y2-y1))/(z2-z1); 

					dx=(xr-x)*Ca-(yr-y)*Sa;
					dy=(xr-x)*Sa+(yr-y)*Ca;

					d2 = dx*dx+dy*dy;
					if(d2 > maxd2){  
						maxd2 = d2;
					}

					s+=((double)abs((int)dx)*(z-oldz));

					nr++;

					oldz=z;
				}
			}

			if (oldz>z1) {
				s/=(oldz-z1);
			}
			else {
				s=0;
			}



			if (s>0 && s<mins) {
				mins=s;
				if (k==0) etaMin=eta;
				else etaMin=-eta;
			}

			
		}
	}


	if(maxd2 > 0 && (z2 -z1 > 0)){
		kruemmung = sqrt(maxd2);
	}else{
		kruemmung = 0;
	}



    diam__ *= 10;		// 1/10 mm

	if(USE_PERCENTAGE_DIAM_RATION){
		log->kr = (int)((100*(kruemmung))/diam__);  // in % ration of the Diameter
	}else{
		log->kr = (int)kruemmung;  // in 1/10 mm/m
	}

	return TRUE;
  }


  /*

		here I can plug in the invented algorithm:)
		returning dataSetKr.WA == ?
  */
  int hasEinseitigerWurzelAn(LOG *log){
		
	return 0;
  }


 int BestimmeIndex(LOG *log,int lv,int lb,int *iv,int *ib){
	int  i;
	BOOL found;

	found=FALSE;
	for (i=0; i<log->nrslices; i++) {
		if (log->slices[i].z>=lv) {
			*iv=i;
			found=TRUE;
			break;
		}
	}
	if (!found) {
		*iv=max(0,log->nrslices-1);
		*ib=max(0,log->nrslices-1);
		return -1;
	}

	found=FALSE;
	for (; i<log->nrslices; i++) {
		if (log->slices[i].z<=lb) {
			*ib=i;
			found=TRUE;
		} else {
			break;
		}
	}
	if (!found) {
		*ib=*iv;
		return -2;
	}

	return 0;


}


 BOOL getKrAngle(LOG *log, int off_start, int off_end, int diam__, double *dWinkel){
 	double x,y,z,xr,yr,x1,y1,z1,x2,y2,z2,oldz,dx,dy;
	int iv,ib,lv,lb;

	int i,k,nr;

	double eta, etaMin, etaCurv, etaFinal, etaTurn;
	double s, s1, mins;
	double Ca, Sa;
	double kruemmung;
	double d2, maxd2, maxDeg = 0;

	const int size__ = 180+1;
	double map[size__][2];
	double val[size__];
	// Mittelpunkt Anfang
	
	int avgPt = 0;
	avgPt = POS_AVG_OFFSET/2;

	lv= off_start-avgPt;
	lb = off_start+avgPt;

	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x1=0;
	y1=0;
	for (i=iv; i<=ib; i++) {
		x1+=log->slices[i].centerx;
		y1+=log->slices[i].centery;
		nr++;
	}
	if (nr>0) {
		x1/=(double)nr;
		y1/=(double)nr;
	}
	else {
		return FALSE;
	}
	z1=(double)(lv+lb)/2.;


	// Mittelpunkt Ende
	lv= off_end-avgPt;
	lb = off_end+avgPt;

	BestimmeIndex(log,lv,lb,&iv,&ib);

	nr=0;
	x2=0;
	y2=0;
	for (i=iv; i<=ib; i++) {
		x2+=log->slices[i].centerx;
		y2+=log->slices[i].centery;
		nr++;
	}
	if (nr>0) {
		x2/=(double)nr;
		y2/=(double)nr;
	}
	else {
		return FALSE;
	}
	z2=(double)(lv+lb)/2;


	log->rx[0]=(int)x1;
	log->ry[0]=(int)y1;
	log->rz[0]=(int)z1;

	log->rx[1]=(int)x2;
	log->ry[1]=(int)y2;
	log->rz[1]=(int)z2;

	etaMin=0;

	mins=1e10; 

	int degOff = 0;
	for (eta=0.001; eta<=1.58; eta+=0.025) {
		for (k=0; k<2; k++) {
			Sa=sin(eta);    
			if (k==1) Sa=-Sa;
			Ca=cos(eta);
			s=0;  
			oldz=z1;   

			nr=0;
			for (i=0; i<log->nrslices; i++) {
				if (log->slices[i].z>z1 && log->slices[i].z<z2) {

					xr=log->slices[i].centerx;
					yr=log->slices[i].centery;
					z=log->slices[i].z;

					x=x1+((z-z1)*(x2-x1))/(z2-z1);
					y=y1+((z-z1)*(y2-y1))/(z2-z1);  

					dx=(xr-x)*Ca-(yr-y)*Sa;
					dy=(xr-x)*Sa+(yr-y)*Ca;

					d2 = dx*dx+dy*dy;
					
					if(d2 > maxd2){
						maxd2 = d2;
					}

					s+=((double)abs((int)dx)*(z-oldz));

					nr++;

					oldz=z;
				}
			}

			if (oldz>z1) {
				s/=(oldz-z1);
			}
			else {
				s=0;
			}

			if (s>0 && s<mins) {
				mins=s;
				if (k==0) etaMin=eta;
				else etaMin=-eta;
			}
		}

		if(maxd2 > 0 && (z2 -z1 > 0)){
		  kruemmung = sqrt(maxd2);
		}else{
			kruemmung = 0;
		}

	}
///////////////
		eta=etaMin;

	Sa=sin(eta);    
	Ca=cos(eta);
	s=0;  
	s1=0;

	oldz=z1;

	for (i=0; i<log->nrslices; i++) {
		if (log->slices[i].z>z1 && log->slices[i].z<z2) {

			xr=log->slices[i].centerx;
			yr=log->slices[i].centery;
			z=log->slices[i].z;

			x=x1+((z-z1)*(x2-x1))/(z2-z1);
			y=y1+((z-z1)*(y2-y1))/(z2-z1);

			dx=(xr-x)*Ca-(yr-y)*Sa;
			dy=(xr-x)*Sa+(yr-y)*Ca;

			if (dy>0) {
				s+=dy*(z-oldz);   
			}
			else if (dy<0) {                
				s1+=(-dy)*(z-oldz);    
			}  
			oldz=z;
		}

	}       


	if (s1>s) {  
		if (etaMin>0) {
			etaMin-=3.1415;
		}
		else {
			etaMin+=3.1415;
		}
	}

	etaCurv=3.1415/2.0-etaMin;				// ermittelter Kruemmungswinkel
	if (etaCurv<-3.1415) etaCurv+=2.0*3.1415;
	else if (etaCurv>3.1415) etaCurv-=2.0*3.1415;

	etaFinal=(double)90*3.1415/180.;		// gewuenschter Kruemmungswinkel
	if (etaFinal<-3.1415) etaFinal+=2.0*3.1415;
	else if (etaFinal>3.1415) etaFinal-=2.0*3.1415;

	*dWinkel = radMath2Grad(etaCurv);



		/*
		degOff = (int)radMath2GradUhr(eta);

		diam__ *= 10;		// 1/10 mm

		if(USE_PERCENTAGE_DIAM_RATION){
		  kruemmung =   (int)((100*(kruemmung))/diam__);  // in 1/10 mm/m
		}else{
			// already set;
		}

		if(degOff %2 == 0){  // in this case we abstract to 180°
			map[degOff][0] = (int)kruemmung;
			map[degOff][1] = degOff;
			val[degOff] = map[degOff][0];
		}
	}

	qsort(val, size__, sizeof(int), compareMin);

	int max_off = size__/2;
	for(int p = size__-1; p > size__/4; p--){
		if(val[p] != 0){
			max_off = p;
			break;
		}
	}
	maxDeg = val[max_off];  // this is now the max
	
	
	// now get the offset where we found the minimum
	for(int x__ = 0; x__ <= size__; x__++){
		if(map[x__][0] == maxDeg){
			log->krAngle = (int)map[x__][1];
			break;
		}
	}

  */
	return TRUE;
  }


  



#if LOAD_DISHAPE3D_MODELS

BOOL ReadMessungFromFile(char *path, LOG *plog){

	char str[512];
	HANDLE hFile;
	DWORD dwRead;
	int i;
	int kennung;



		
	hFile = CreateFile((LPCTSTR)path,GENERIC_READ ,FILE_SHARE_READ,
						   NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);

	if(hFile == INVALID_HANDLE_VALUE){
		MSG(str);
		return(FALSE);
	}

	ReadFile(hFile,(LPVOID)&kennung,sizeof(kennung),&dwRead,NULL);
	if (kennung==0x6969) {

		ReadFile(hFile,(LPVOID)&plog->la,sizeof(plog->la),&dwRead,NULL);
		ReadFile(hFile,(LPVOID)&plog->nrslices,sizeof(plog->nrslices),&dwRead,NULL);

		if (plog->nrslices<0) plog->nrslices=0;
		if (plog->nrslices>MAX_SLICE) plog->nrslices=0;

		for (i=0; i<plog->nrslices; i++) {
			ReadFile(hFile,(LPVOID)&plog->slices[i],sizeof(plog->slices[i]),&dwRead,NULL);

		}
	} else {
		SetFilePointer(hFile,0,NULL,FILE_BEGIN);
		ReadFile(hFile,(LPVOID)plog,sizeof(LOG),&dwRead,NULL);
	}

	
	plog->loffset=0;
	CloseHandle(hFile);

	
	return TRUE;

}

BOOL ReadMessung(char *path, LOG *plog)
{
	char lpszPath[_MAX_PATH];
	char str[512];
	HANDLE hFile;
	DWORD dwRead;

	int i;
	int kennung;


	return ReadMessungFromFile(path, plog);

}

#endif


