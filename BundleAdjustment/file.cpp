#include "stdafx.h"
#include "file.h"

//////////////////////////////
// 特徴点入力ファイル読み込み
//
void readfilePoints(char * filename, Mat& points, Mat& imageX, Mat& imageY)
{
	///---ファイル読み込み---
	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "r");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	char str[SHRT_MAX];
	char * delim = " ";
	char * ctx;
	char * token;
	int i=0, j=0;
	///---1行づつ読み込んでいく---
	while( fgets(str, LONG_MAX, fp) != NULL ){
		token = strtok_s( str, delim, &ctx); //半角スペースで分割
		while( token != NULL ){
			//分割した文字列をatofしてcorrpointに代入していく
			if(j==0) points.at<float>(i,0) = (float)atof(token);
			else if(j==1) points.at<float>(i,1) = (float)atof(token);
			else if(j==2) points.at<float>(i,2) = (float)atof(token);
			else if(j%2==1 && j<(3+FRAMETOTAL*2)) imageX.at<float>(i,(j-3)/2) = (float)atof(token);
			else if(j%2==0 && j<(3+FRAMETOTAL*2)) imageY.at<float>(i,(j-4)/2) = (float)atof(token);
			//else //何もしない

			token = strtok_s(NULL, delim, &ctx);  //2回目以降
			j++;
		}
		i++;
		j=0;
	}
	fclose( fp );
}

////////////////////////////
// CorrPoint構造体初期化
//
void initCorrPoint(CorrPoint * corrpoint)
{
	for(int i=0;i<KPtotal;i++)
	{
		corrpoint[i].x0 = (float *)malloc(sizeof(float)*FRAMETOTAL);
		corrpoint[i].y0 = (float *)malloc(sizeof(float)*FRAMETOTAL);
	}
}

////////////////////////////
// CorrPoint構造体free
//
void freeCorrPoint(CorrPoint * corrpoint)
{
	for(int i=0;i<KPtotal;i++)
	{
		free(corrpoint[i].x0);
		free(corrpoint[i].y0);
	}
}

/////////////////////////////
//特徴点入力ファイル読み込み
//
void readfileCorr(char * filename, CorrPoint * corrpoint)
{
	///---ファイル読み込み---
	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "r");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	char str[SHRT_MAX];
	char * delim = " ";
	char * ctx;
	char * token;
	int i=0, j=0;
	///---1行づつ読み込んでいく---
	while( fgets(str, LONG_MAX, fp) != NULL ){
		token = strtok_s( str, delim, &ctx); //半角スペースで分割
		while( token != NULL ){
			//分割した文字列をatofしてcorrpointに代入していく
			if(j==0) corrpoint[i].x = (float)atof(token);
			else if(j==1) corrpoint[i].y = (float)atof(token);
			else if(j==2) corrpoint[i].z = (float)atof(token);
			else if(j%2==1 && j<(3+FRAMETOTAL*2)) corrpoint[i].x0[(j-3)/2] = (float)atof(token);
			else if(j%2==0 && j<(3+FRAMETOTAL*2)) corrpoint[i].y0[(j-4)/2] = (float)atof(token);
			//else //何もしない

			token = strtok_s(NULL, delim, &ctx);  //2回目以降
			j++;
		}
		i++;
		j=0;
	}
	fclose( fp );
}

///////////////////////////////
// ファイル読み捨てて行数取得
//
int readfileLine(char * filename)
{
	int total = 0;

	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "r");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	char s[SHRT_MAX];
	while( fgets( s, SHRT_MAX, fp ) != NULL ){
		total++;
	}

	return total;
}

/////////////////////////////
//並進・回転ファイル読み込み
//
void readfileRT(char * filename, Mat * rt)
{
	// Open File Storage
	cv::FileStorage cvfs(filename,CV_STORAGE_READ);

	cv::FileNode node(cvfs.fs, NULL);	// Get Top Node
	cv::FileNode fn = node["mat_rt"]; //ノードを指定

	for(int i=0;i<FRAMETOTAL;i++){
		cv::read(fn[i], rt[i]);
	}
}

/////////////////////////////////
//カメラ内部行列ファイル読み込み
//
void readfileCamera(char * filename, Mat * cameraA)
{
	// Open File Storage
	cv::FileStorage cvfs(filename,CV_STORAGE_READ);

	cv::FileNode node(cvfs.fs, NULL);	// Get Top Node
	cv::FileNode fn = node["intrinsic"]; //ノードを指定

	cv::read(fn, *cameraA);
}

///////////////////////////
// paraAファイル書き込み
//
void writefile(char * filename, float * paraA, WriteType type , float min, float max, Mat * rt, CorrPoint * corrpoint)
{
	if(type==0)//XYZ
	{
		///---ファイル読み込み---
		FILE * fp;
		errno_t error;
		error = fopen_s(&fp, filename, "w");
		if(error != 0){
			cout << "ファイルが開けません " << filename << endl;
			exit(1);
		}

		//色つきかどうか
		if(min == -1. && max == -1.){
			for(int i=0;i<KPtotal;i++){
				fprintf(fp,"%f ",paraA[i*3+1]);
				fprintf(fp,"%f ",paraA[i*3+2]);
				fprintf(fp,"%f",paraA[i*3+3]);
				fprintf(fp,"\n");
			}
		}else{
			int white = 255*256*256 + 255*256 + 255;
			int arrowNum = 20;
			//カメラ位置表示するかどうか
			if(rt != NULL){
				fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", KPtotal+FRAMETOTAL*arrowNum, KPtotal+FRAMETOTAL*arrowNum);
				for(int i=0;i<FRAMETOTAL;i++){
					/*
					Mat camera = Mat::eye(4,4,CV_32F);
					Mat temp = Mat::eye(4,4,CV_32F);
					for(int j=0;j<i+1;j++){
					Mat r = rt[j].operator()(Range(0,3),Range(0,3));
					Mat t = rt[j].operator()(Range(3,4),Range(0,3));
					temp.at<float>(0,0) = r.at<float>(0,0);
					temp.at<float>(0,1) = r.at<float>(0,1);
					temp.at<float>(0,2) = r.at<float>(0,2);
					temp.at<float>(1,0) = r.at<float>(1,0);
					temp.at<float>(1,1) = r.at<float>(1,1);
					temp.at<float>(1,2) = r.at<float>(1,2);
					temp.at<float>(2,0) = r.at<float>(2,0);
					temp.at<float>(2,1) = r.at<float>(2,1);
					temp.at<float>(2,2) = r.at<float>(2,2);
					temp.operator()(Range(0,3),Range(3,4)) = t.t();
					camera = camera * temp;
					}
					//cout << camera << endl;

					for(int j=0;j<=arrowNum;j++){
					fprintf(fp,"%f ",camera.at<float>(3,0)+(float)(arrowNum-j)/arrowNum*camera.at<float>(2,0)*10);
					fprintf(fp,"%f ",camera.at<float>(3,1)+(float)(arrowNum-j)/arrowNum*camera.at<float>(2,1)*10);
					fprintf(fp,"%f ",camera.at<float>(3,2)+(float)(arrowNum-j)/arrowNum*camera.at<float>(2,2)*10);
					fprintf(fp,"%d",white);
					fprintf(fp,"\n");
					}
					*/
					for(int j=0;j<=arrowNum;j++){
						fprintf(fp,"%f ",rt[i].at<float>(3,0)+(float)(arrowNum-j)/arrowNum*rt[i].at<float>(0,2)*10);
						fprintf(fp,"%f ",rt[i].at<float>(3,1)+(float)(arrowNum-j)/arrowNum*rt[i].at<float>(1,2)*10);
						fprintf(fp,"%f ",rt[i].at<float>(3,2)+(float)(arrowNum-j)/arrowNum*rt[i].at<float>(2,2)*10);
						fprintf(fp,"%d",white-40*i);
						fprintf(fp,"\n");
					}
				}
			}else{
				fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", KPtotal, KPtotal);
			}
			int red = 255*256*256;
			if(corrpoint == NULL){
				for(int i=0;i<KPtotal;i++){
					fprintf(fp,"%f ",paraA[i*3+1]);
					fprintf(fp,"%f ",paraA[i*3+2]);
					fprintf(fp,"%f ",paraA[i*3+3]);
					fprintf(fp,"%d",(int)(red + 256*(int)(255*(1.-(paraA[i*3+3]-min)/(max-min)))));
					fprintf(fp,"\n");
				}
			}else{
				for(int i=0;i<KPtotal;i++){
					fprintf(fp,"%f ",paraA[i*3+1]);
					fprintf(fp,"%f ",paraA[i*3+2]);
					fprintf(fp,"%f ",paraA[i*3+3]);
					int frame=0;
					for(int j=0;j<FRAMETOTAL;j++){
						if(corrpoint[i].x0[j] != -999.){
							frame = j;
							break;
						}
					}
					fprintf(fp,"%d",(int)(red + 256*(int)(255*((float)frame/(FRAMETOTAL-1)))));
					fprintf(fp,"\n");
				}
			}
		}
		fclose(fp);
	}
}

///////////////////
//Mat書き出し
//
void writefileMat(char * filename, Mat * mat)
{
	// Open File Storage
	cv::FileStorage	cvfs(filename,CV_STORAGE_WRITE);
	cv::WriteStructContext ws(cvfs, "mat", CV_NODE_SEQ);	// create node
	for(int i=0; i<FRAMETOTAL; i++){
		cv::write(cvfs,"",mat[i]);
	}
}

///////////////////
//Mat書き出し
//
void writefileDistort(char * filename, Mat Distort)
{
	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	fprintf(fp,"K1 K2 K3 P1 P2\n\n");
	for(int i=0;i<FRAMETOTAL;i++){
		float K1 = Distort.at<float>(0,i);
		float K2 = Distort.at<float>(1,i);
		float P1 = Distort.at<float>(2,i);
		float P2 = Distort.at<float>(3,i);
		float K3 = Distort.at<float>(4,i);
		fprintf(fp,"フレーム%d\n",i);
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\n\n",K1,K2,K3,P1,P2);
	}
	fclose(fp);
}


///////////////////////////
// 3次元点群書き込み
//
void writefilePoints(char * filename, Mat& points, int pointTotal)
{
	///---ファイル読み込み---
	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}

	for(int i=0;i<pointTotal;i++){
		fprintf(fp,"%f,",points.at<float>(i,0));
		fprintf(fp,"%f,",points.at<float>(i,1));
		fprintf(fp,"%f",points.at<float>(i,2));
		fprintf(fp,"\n");
	}
	// */

	fclose(fp);
}

///////////////////////////
// 3次元点群書き込み
//
void writefilePoints(char * filename, Mat& points, int pointTotal, int frameTotal, float min, float max, Mat * R, Mat * T)
{
	///---ファイル読み込み---
	FILE * fp;
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	
	int red = 255*256*256;
	int green = 255*256*256 + 255*256;
	int white = 255*256*256 + 255*256 + 255;
	int arrowNum = 20;

	/*
	if(R != NULL){
		///カメラ位置表示する
		fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", pointTotal+frameTotal*(arrowNum+1), pointTotal+frameTotal*(arrowNum+1));
		for(int i=0;i<frameTotal;i++){
			for(int j=0;j<=arrowNum;j++){
				fprintf(fp,"%f ",T[i].at<float>(0,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(0,2)*5);
				fprintf(fp,"%f ",T[i].at<float>(1,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(1,2)*5);
				fprintf(fp,"%f ",T[i].at<float>(2,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(2,2)*5);
				fprintf(fp,"%d",white-40*i);
				fprintf(fp,"\n");
			}
		}
	}else{
		///カメラ位置表示しない
		fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", pointTotal, pointTotal);
	}

	for(int i=0;i<pointTotal;i++){
		fprintf(fp,"%f ",points.at<float>(i,0));
		fprintf(fp,"%f ",points.at<float>(i,1));
		fprintf(fp,"%f ",points.at<float>(i,2));
		fprintf(fp,"%d",(int)(red + 256*(int)(255*(1.-(points.at<float>(i,2)-min)/(max-min)))));
		fprintf(fp,"\n");
	}
	// */
	
	///*
	if(R != NULL){
		///カメラ位置表示する
		fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", pointTotal+frameTotal*(arrowNum+1), pointTotal+frameTotal*(arrowNum+1));
		for(int i=0;i<frameTotal;i++){
			for(int j=0;j<=arrowNum;j++){
				fprintf(fp,"%f ",T[i].at<float>(0,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(0,2)*5);
				fprintf(fp,"%f ",T[i].at<float>(1,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(1,2)*5);
				fprintf(fp,"%f ",T[i].at<float>(2,0)+(float)(arrowNum-j)/arrowNum*R[i].at<float>(2,2)*5);
				fprintf(fp,"%d",white);
				fprintf(fp,"\n");
			}
		}
	}else{
		///カメラ位置表示しない
		fprintf(fp,"# .PCD v.7 - Point Cloud Data file format\nVERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS %d\nDATA ascii\n", pointTotal, pointTotal);
	}

	for(int i=0;i<pointTotal;i++){
		fprintf(fp,"%f ",points.at<float>(i,0));
		fprintf(fp,"%f ",points.at<float>(i,1));
		fprintf(fp,"%f ",points.at<float>(i,2));
		fprintf(fp,"%d",white);
		fprintf(fp,"\n");
	}
	// */

	fclose(fp);
}