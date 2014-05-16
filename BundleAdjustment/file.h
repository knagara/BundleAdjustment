#ifndef FILE_H_INCLUDED
#define FILE_H_INCLUDED

#include "stdafx.h"
#include "opencv.h" //OpenCVライブラリ

extern const int FRAMETOTAL,WID,HEI;
extern int KPtotal;

///---入力データを入れるための構造体CorrPoint---
//CorrPointとは「Corresponding Point = 対応点」のことですが特徴点と同じだと思ってください。複数フレーム間で対応している点なので対応点と呼んでいます
typedef struct{
	float x;
	float y;
	float z;
	float * x0;
	float * y0;
}CorrPoint;

///---write関数で使用する列挙型---
enum WriteType{
	XYZ,
	R,
	T,
	FXFYCXCY
};

///---関数プロトタイプ宣言---
void readfilePoints(char * filename, Mat& points, Mat& imageX, Mat& imageY);
void initCorrPoint(CorrPoint * corrpoint);
void freeCorrPoint(CorrPoint * corrpoint);
void readfileCorr(char * filename, CorrPoint * corrpoint);
int readfileLine(char * filename);
void readfileRT(char * filename, Mat * rt);
void readfileCamera(char * filename, Mat * cameraA);
void writefile(char * filename, float * paraA, WriteType type, float min=-1., float max=-1., Mat * rt = NULL, CorrPoint * corrpoint = NULL);
void writefileMat(char * filename, Mat * rt);
void writefileDistort(char * filename, Mat Distort);
void writefilePoints(char * filename, Mat& points, int pointTotal);
void writefilePoints(char * filename, Mat& points, int pointTotal, int frameTotal, float min, float max, Mat * R = NULL, Mat * T = NULL);

#endif