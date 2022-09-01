// completeGraph.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "completeGraph.h"
#include "cv.h"
#include "highgui.h"
#include "cvaux.h"
#include "cxcore.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstddef>
#include <math.h>

extern "C"
{
	#include "imgfeatures.h"
    #include "utils.h"
}


using namespace std;
using namespace cv;

typedef struct  _tagTiePoint
{
	long ptID;
	float xme;
	float yme;
	float xother;
	float yother;
}TIEPOINT;

typedef struct _tagSurpotFactor
{
	double factor;
	TIEPOINT tiePt;
}SF;
#define E 2.718281828459
void getTriangleDescriptor( double *tp1, double *tp2, double *tp3, double * descriptor);
double get3DVectorLength(double *vec);
int cmpSfs(const void *a, const void *b);
vector<SF> grossErrorDetectionbyTriangleConstraint( vector<feature> feats);
int compare(const void *a,const void *b)
{
	if(*(int*)b-*(int*)a>0)
		return 1;
	else return -1;
	return 0;
}


int getGridWidthAndHeight(const char *path ,int *fxlStart , int *fylStart ,int *fxlEnd , int *fylEnd)
{
	ifstream inf(path,ios::in);

	if (!inf.is_open())
	{
		cout<<"file open failed!"<<endl;
		return -1;
	}	
	int numL=0 , numR=0;
	
	int t;
	inf>>numL>>numR>>t>>t;
	int *fxl=new int[numL] ; int *fyl=new int[numL];
	int *fxr=new int[numR] ; int *fyr=new int[numR];
	double flag;

	double xl,yl,xr,yr;
	for (int i=0 ; i<numL ; ++i)
	{
		inf>>t>>xl>>yl>>xr>>yr>>flag>>flag>>t; 
		fxl[i] = int(xl * 118.110 + 0.5);
		fyl[i] = int(yl * 118.110 + 0.5);
		fxr[i] = int(xr * 118.110 + 0.5);
		fyr[i] = int(yr * 118.110 + 0.5);
	}
	inf.close();

	FILE *fp = fopen("TiePts.dpt","w");
	fprintf(fp,"%d %d 0\n",numL,numR);
	for (int i=0;i<numL;++i)
	{
		fprintf(fp,"%d %d %d %d %d 0.0 0.0 1\n",i,fxl[i],400-fyl[i],fxr[i],400-fyr[i]);
	}
	fclose(fp);

	qsort(fxl,numL,sizeof(fxl[0]),compare);
	qsort(fyl,numL,sizeof(fyl[0]),compare);

	*fxlStart = fxl[numL - 1];
	*fylStart = fyl[numL - 1];
	*fxlEnd = fxl[0];
	*fylEnd = fyl[0];

	delete []fxl; delete []fyl; delete []fxr; delete []fyr;
	 
	return 1;
}
int getGrid( const char *path )
{
	const double gridSz = 400.0;

	int maxx,maxy,minx,miny;
	getGridWidthAndHeight(path ,&minx ,&miny ,&maxx ,&maxy );
	int numLx = floor((maxx - minx)/gridSz + 0.5);
	int numLy = floor((maxy - miny)/gridSz + 0.5);

	ifstream inf(path,ios::in);

	if (!inf.is_open())
	{
		cout<<"file open failed!"<<endl;
		return -1;
	}	
	int numL=0 , numR=0;
	int t;
	inf>>numL>>numR>>t>>t;
	int *fxl=new int[numL] ; int *fyl=new int[numL];
	int *fxr=new int[numR] ; int *fyr=new int[numR];
	double flag;
	
	double xl,yl,xr,yr;
	for (int i=0 ; i<numL ; ++i)
	{
		inf>>i>>xl>>yl>>xr>>yr>>flag>>flag>>t; 
		fxl[i] = xl * 118.110;
		fyl[i] = yl * 118.110;
		fxr[i] = xr * 118.110;
		fyr[i] = yr * 118.110;
	}
	inf.close();

	
	
	for(int m=0 ; m<numLx ; m++)
	{
		double startx = m*gridSz;
		for(int n=0 ; n<numLy ; n++)
		{
			double starty = n*gridSz;
			vector<feature> feats;
			for(int i=0 ; i<numL ; i++)
			{
				if(fxl[i]>=startx && fxl[i]<startx + gridSz && fyl[i] >= starty && fyl[i]<starty + gridSz)
				{
					struct feature feat;
					feat.x=fxl[i];
					feat.y=fyl[i];
					feat.xother=fxr[i];
					feat.yother=fyr[i];
					feats.push_back(feat);
				}  
			}
			int nValidPts=0;
			int nCounter = 0;
			if (feats.size() <= 3 )
				continue;
			do 
			{
				//vector<SF> sfs=grossErrorDetectionbyTriangleConstraint(feats , feats.size() , nValidPts);
				//nCounter ++;
			} while (nCounter <= 50);
		}
	}
				
	
	//double a=feats[3]->x;
	return 1;
}
vector<SF> grossErrorDetectionbyTriangleConstraint( vector<feature> feats)
{
	int n = feats.size();
	vector<SF> SFs(n);//attributes of graph nodes
	TIEPOINT tp1,tp2,tp3;
	for (int i=0;i<n;++i)
	{
		double factor = 0;
		int count = 0;
		tp1.xme = feats[i].x;
		tp1.yme = feats[i].y;
		tp1.xother = feats[i].xother;
		tp1.yother = feats[i].yother;
		tp1.ptID = feats[i].d;
		for (int j = 0; j<n;++j)
		{
			tp2.xme = feats[j].x;
			tp2.yme = feats[j].y;
			tp2.xother = feats[j].xother;
			tp2.yother = feats[j].yother;
			for (int k=0;k<n;++k)
			{
				tp3.xme = feats[k].x;
				tp3.yme = feats[k].y;
				tp3.xother = feats[k].xother;
				tp3.yother = feats[k].yother;
				if (i != j && j != k && i != k)
				{
					double dx1 = tp2.xme - tp1.xme, dy1 = tp2.yme - tp1.yme;
					double dx2 = tp3.xme - tp2.xme, dy2 = tp3.yme - tp2.yme;
					double dx3 = tp3.xme - tp1.xme, dy3 = tp3.yme - tp1.yme;
					if(fabs(dx2*dy1 - dy2*dx1) <= 0.1)//approximately parallel
						continue;
					if (sqrt(dx1*dx1+dy1*dy1) < 3.0 || sqrt(dx2*dx2+dy2*dy2) < 3.0 || sqrt(dx3*dx3+dy3*dy3) < 3.0)//features are two close
						continue;

					dx1 = tp2.xother - tp1.xother; dy1 = tp2.yother - tp1.yother;
					dx2 = tp3.xother - tp1.xother; dy2 = tp3.yother - tp1.yother;
					dx3 = tp3.xother - tp2.xother; dy3 = tp3.yother - tp2.yother;
					if(fabs(dx2*dy1 - dy2*dx1) <= 0.1)
						continue;
					if (sqrt(dx1*dx1+dy1*dy1) < 3.0 || sqrt(dx2*dx2+dy2*dy2) < 3.0 || sqrt(dx3*dx3+dy3*dy3) < 3.0)
						continue;

					double descriptor1[3] = {0}, descriptor2[3] = {0};
					double lpt1[2] = {tp1.xme,tp1.yme}, lpt2[2] = {tp2.xme,tp2.yme}, lpt3[2] = {tp3.xme,tp3.yme};
					double rpt1[2] = {tp1.xother,tp1.yother}, rpt2[2] = {tp2.xother,tp2.yother}, rpt3[2] = {tp3.xother,tp3.yother};
					getTriangleDescriptor(lpt1,lpt2,lpt3,descriptor1);
					getTriangleDescriptor(rpt1,rpt2,rpt3,descriptor2);
					double vec[3] = {descriptor1[0]-descriptor2[0],descriptor1[1]-descriptor2[1],descriptor1[2]-descriptor2[2]};
					double length = get3DVectorLength(vec);

					factor += 1.0/pow(E,length);
					count ++;
				}
			}
		}
		SFs[i].tiePt = tp1;
		SFs[i].factor = factor/count;
	}
	qsort(&(SFs[0]),SFs.size(),sizeof(SF),cmpSfs);

	return SFs;
}

inline void getTriangleDescriptor( double *tp1, double *tp2, double *tp3, double * descriptor)
{
	double vec12[2] = {tp2[0]-tp1[0],tp2[1]-tp1[1]};
	double vec13[2] = {tp3[0]-tp1[0],tp3[1]-tp1[1]};
	double vec21[2] = {-vec12[0],-vec12[1]};
	double vec23[2] = {tp3[0]-tp2[0],tp3[1]-tp2[1]};
	double vec31[2] = {-vec13[0],-vec13[1]};
	double vec32[2] = {-vec23[0],-vec23[1]};

	double length12 = sqrt(vec12[0]*vec12[0]+vec12[1]*vec12[1]);
	double length13 = sqrt(vec13[0]*vec13[0]+vec13[1]*vec13[1]);
	double length23 = sqrt(vec23[0]*vec23[0]+vec23[1]*vec23[1]);

	descriptor[0] = (vec12[0]*vec13[0]+vec12[1]*vec13[1])/(length12*length13);
	descriptor[1] = (vec21[0]*vec23[0]+vec21[1]*vec23[1])/(length23*length12);
	descriptor[2] = (vec31[0]*vec32[0]+vec31[1]*vec32[1])/(length13*length23);
}

inline double get3DVectorLength( double *vec )
{
	return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}
int cmpSfs( const void *a, const void *b )
{
	SF *sfa = (SF*)a, *sfb = (SF*)b;
	if ((*sfa).factor - (*sfb).factor < 0)
		return 1;
	else if ((*sfa).factor - (*sfb).factor > 0)
		return -1;
	else
		return 0;
}

void deleteElementbyID(vector<feature> & features, int targetID)
{
	int n = features.size();
	for (int i=0;i<n;++i)
	{
		if (features[i].d == targetID)
		{
			features.erase(features.begin()+i);
			break;
		}
	}
}

struct feature getElementbyID(vector<feature> &features, int targetID)
{
	int n = features.size();
	for (int i=0;i<n;++i)
	{
		if (features[i].d == targetID)
			return features[i];
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	const char *filePath="TiePts.dpt";
	vector<Point2f> ptsL,ptsR;
	FILE *fp = fopen(filePath,"r");
	int id,xl,yl,xr,yr,flag;
	float f;
	vector<feature> ori_feats;
	while(8 == fscanf(fp,"%d%d%d%d%d%f%f%d",&id,&xl,&yl,&xr,&yr,&f,&f,&flag))
	{
		ptsL.push_back(Point2f(xl,yl));
		ptsR.push_back(Point2f(xr,yr));
		struct feature feat;
		feat.d = id;
		feat.x=xl;
		feat.y=yl;
		feat.xother=xr;
		feat.yother=yr;
		ori_feats.push_back(feat);
	}

	RNG rng;
	int *rng_arr = new int[10000];
	rng.fill(Mat(10000,1,CV_32SC1,rng_arr),RNG::UNIFORM,0,400);

	int nGrossErrors = 100;
	for (int i=0;i<nGrossErrors;++i)
	{
		struct feature feat;
		feat.d = 1000+i;
		feat.x = rng_arr[i*4+0];
		feat.y = rng_arr[i*4+1];
		feat.xother = rng_arr[i*4+2];
		feat.yother = rng_arr[i*4+3];
		ori_feats.push_back(feat);
	}

	vector<feature> detected_feats1,detected_feats2;
	double c1,e1,c2,e2;

	{//Region: error detection by the triangle constraint.
		vector<feature> t_feats = ori_feats;
		vector<SF> sfs = grossErrorDetectionbyTriangleConstraint(t_feats);
		double last_min_Attribut = sfs[sfs.size()-1].factor;;
		double ds = 0;
		do 
		{
			int id = sfs[sfs.size()-1].tiePt.ptID;
			detected_feats1.push_back(getElementbyID(t_feats,id));
			deleteElementbyID(t_feats,id);
			sfs = grossErrorDetectionbyTriangleConstraint(t_feats);
			double now_min_Attribute = sfs[sfs.size()-1].factor;
			ds = abs(last_min_Attribut-now_min_Attribute);
			last_min_Attribut = now_min_Attribute;
		} while (ds > 0.01 || last_min_Attribut < 0.77);

		int nc = 0, ne = 0;
		for (int i=0;i<detected_feats1.size();++i)
		{
			if (detected_feats1[i].d < 1000)
				nc++;
			else
				ne++;
		}
		c1 = nc/30.0;
		e1 = 1.0*ne/nGrossErrors;
	}

	{//Region: error detection by the homograhpy constraint.
		vector<Point2f> t_ptsL,t_pstR;
		for (int i=0;i<ori_feats.size();++i)
		{
			double t_xl = ori_feats[i].x;
			double t_yl = ori_feats[i].y;
			double t_xr = ori_feats[i].xother;
			double t_yr = ori_feats[i].yother;
			t_ptsL.push_back(Point2f(t_xl,t_yl));
			t_pstR.push_back(Point2f(t_xr,t_yr));
		}
		uchar *mask_arr = new uchar[ori_feats.size()];
		findHomography(t_ptsL,t_pstR,CV_RANSAC,3.0,Mat(ori_feats.size(),1,CV_8UC1,mask_arr));
		int fm = 0;
		int nc = 0,  ne = 0;
		for (int i=0;i<t_ptsL.size();++i)
		{
			if (mask_arr[i] == 0 && ori_feats[i].d >= 1000)
				ne++;
			if (mask_arr[i] == 0 && ori_feats[i].d < 1000)
				nc ++;
		}
		c2 = nc/30.0;
		e2 = 1.0*ne/nGrossErrors;
	}
	


	 //getGrid(filePath);
	return 0;
}


