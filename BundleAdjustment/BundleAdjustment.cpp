// BundleAdjustment.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include "opencv.h" //OpenCVライブラリ
#include "file.h" //ファイルを扱う関数群
#include "bundle.h" //バンドル調整


////////////////////////
// グローバル変数宣言
//
selfParam selfFlag = SELF_CALIB_ON;
selfAlterParam selfAlterFlag = ALTERNATE_OFF;
fixParam fixFlag = FIX_7_AUTO;
fcParam fcFlag = FC_FIX;
float c = 0.0001F, epsilon = 0.01F;
const int FRAMETOTAL = 5; //フレーム数
const int BNDL = FRAMETOTAL - 1; //バンドル幅（フレーム数 - 1）
const int WID = 128, HEI = 128; //画像サイズ
const char filedir[] = "data"; //データを置くディレクトリ名
const char currentdir[] = "sample"; //その中で使用するディレクトリ名 - 最終的に ./filedir/currentdir/hoge.txt のような形になる

const char corrfilename[] = "allbundlepoints.csv"; //特徴点入力データ
const char rtfilename[] = "rt.xml"; //並進・回転ベクトル入力データ
const char camerafilename[] = "camera.xml"; //カメラ内部行列入力データ
const char corroutfilename[] = "result_xyz"; //特徴点出力データ
const char rtoutfilename[] = "result_rt.xml"; //並進・回転ベクトル出力
const char Routfilename[] = "result_R.xml"; //並進・回転ベクトル出力
const char Toutfilename[] = "result_T.xml"; //並進・回転ベクトル出力
const char cameraoutfilename[] = "result_camera.xml"; //カメラ内部行列出力
const char distortoutfilename[] = "result_distort.txt"; //カメラ歪み系数出力
int KPtotal; //特徴点の実際の数「KP = Key Point = 特徴点」
int KP2Dtotal; //特徴点がフレーム毎に写っている数を合計したもの
float fx, fy, cx, cy; //キャリブレーションで求まったカメラ内部パラメータ
//double PI = 3.14159265358979;
//int fps = 15; 

////クイックソート比較関数
int comp( const void *c1, const void *c2 );


//////////////////
// main関数
//
int _tmain(int argc, _TCHAR* argv[])
{
	cout << "-------------------------" << endl;
	cout << "    Bundle Adjustment    " << endl;
	cout << "-------------------------" << endl;
	clock_t start_time_total,end_time_total;
	start_time_total = clock(); //実行時間計測


	///---特徴点入力データ読み込み---
	char filename[255];
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,corrfilename);
	//cout << "特徴点入力データ" << filename << " 読み込みます" << endl;
	//ファイルの行数(=特徴点の数)取得してKPtotalに代入
	KPtotal = readfileLine(filename); //readfileLine関数 - file.cpp参照
	//cout << "完了　特徴点の数 = " << KPtotal << "個" << endl << endl;
	//特徴点入れるMat 画像座標入れるMat
	Mat points = Mat::zeros(KPtotal,3,CV_32F);
	Mat imageX = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F);
	Mat imageY = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F);
	//ファイル読み込んでMatに代入
	readfilePoints(filename, points, imageX, imageY);
	
	//構造体CorrPointインスタンス化 - 構造体の宣言はfile.hに書いてあります
	CorrPoint * corrpoint;
	corrpoint = (CorrPoint *)malloc(sizeof(CorrPoint) * KPtotal);
	initCorrPoint(corrpoint); //corrpointの初期化 - file.cpp参照
	//ファイル読み込んでcorrpointに代入
	readfileCorr(filename, corrpoint); //readfileCorr関数 - file.cpp参照
	


	///---カメラ並進・回転ベクトル読み込み---
	/////ICPnormalMultiDataでrt.xmlを出力した場合，RT行列の形に注意すること！通常の4*4RT行列の形とは違う。
	/////一番下の行は[0,0,0,1]ではなく[tx,ty,tz,1]になっている。
	Mat rt[FRAMETOTAL]; //RT4*4行列入れるMat作成
	for(int i=0;i<FRAMETOTAL;i++) rt[i] = Mat_<float>(4,4); //4*4のメモリ確保
	//rt.xmlファイル読み込んでMatに代入
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,rtfilename);
	//cout << "並進・回転データ" << filename << " 読み込みます" << endl;
	readfileRT(filename, rt); //readfileRT関数 - file.cpp参照
	//cout << "完了" << endl << endl;


	///---カメラ内部行列(fx,fy,cx,cy)読み込み---
	////OpenCVのカメラキャリブレーションで求まる(cx,cy)は画像左上中心ですが，
	////このプログラムでは，(cx,cy)は画像中心からのズレを表します。
	////少し先で，中心からのズレに変換しています
	////ここではデータを読み込むだけです
	Mat cameraA(3,3,CV_32F); //カメラ内部行列 3*3行列 [fx, 0, cx; 0, fy, cy; 0, 0, 1]
	//camera.xmlファイル読み込んでMatに代入
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,camerafilename);
	//cout << "カメラ内部行列データ" << filename << " 読み込みます" << endl;
	readfileCamera(filename, &cameraA); //readfileCamera関数 - file.cpp参照
	//cout << "完了" << endl << endl;


	///---Bundleに渡すためにデータ整形---
	Mat R[FRAMETOTAL]; //R3*3行列入れるMat作成
	Mat T[FRAMETOTAL]; //R3*1行列入れるMat作成
	Mat K[FRAMETOTAL]; //K3*3行列入れるMat作成
	for(int i=0;i<FRAMETOTAL;i++){
		R[i] = Mat_<float>(3,3); //3*3のメモリ確保
		R[i] = rt[i](Range(0,3),Range(0,3));
		T[i] = Mat_<float>(3,1); //3*1のメモリ確保
		T[i] = rt[i](Range(3,4),Range(0,3)).t();
		K[i] = Mat_<float>(3,3); //3*3のメモリ確保
		cameraA.copyTo(K[i]); //深いコピー
	}


	///---バンドル調整---
	//cout << endl;
	//cout << "------------------------------" << endl;
	//cout << "   バンドル調整を開始します   " << endl;
	//cout << "------------------------------" << endl;
	clock_t start_time_bundle,end_time_bundle;
	start_time_bundle = clock(); //実行時間計測


	///---クラスBundleのインスタンス化---
	Bundle bundle(points, imageX, imageY, R, T, K, WID, HEI);
	///---条件設定---
	//bundle.set_imageNaN(-1.); //カメラに投影されてないときの画像座標
	bundle.set_c(c);     //レーベンバーグ・マーカート法のc
	bundle.set_epsilon(epsilon);    //終了条件のイプシロン
	///---バンドル調整start---
	bundle.start(selfFlag, fixFlag, fcFlag, selfAlterFlag);
	

	end_time_bundle = clock(); //実行時間計測終了
	cout << "バンドル調整実行時間 = " << (float)(end_time_bundle - start_time_bundle)/CLOCKS_PER_SEC << "秒" << endl << endl;

	
	///---点群(x,y,z)の出力---
	///(x,y,z)のcsvファイル
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.csv", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows); //writefile関数 - file.cpp参照
	///PointCloudLibraryの点群データ形式.pcd
	//zの値によって　近⇔遠　緑⇔赤　と色指定する
	//zの最大値と最小値を求める
	double minVal, maxVal;
	Point minLoc, maxLoc;
	minMaxLoc(points(Range(0,points.rows),Range(2,3)), &minVal, &maxVal, &minLoc, &maxLoc);
	//ファイル書き込み
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.pcd", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows, imageX.cols, (float)minVal, (float)maxVal, R, T); //writefile関数 - file.cpp参照
	

	///---回転Rの出力---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Routfilename);
	writefileMat(filename, R); //writefileMat関数 - file.cpp参照
	

	///---並進Tの出力---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Toutfilename);
	writefileMat(filename, T); //writefileMat関数 - file.cpp参照

	
	///---カメラ内部行列(fx,fy,cx,cy)の出力---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,cameraoutfilename);
	writefileMat(filename, K); //writefileMat関数 - file.cpp参照

	
	///---カメラ歪み係数(K1,K2,P1,P2,K3)の出力---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,distortoutfilename);
	writefileDistort(filename, bundle.Distort); //writefileMat関数 - file.cpp参照

	
	///---再投影誤差vectorE書き込み---
	FILE * fp;
	sprintf_s(filename, _countof(filename), "%s/%s/result_E.csv",filedir,currentdir);
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "ファイルが開けません " << filename << endl;
		exit(1);
	}
	for(int i=0;i<bundle.vectorE.size();i++){
		fprintf(fp,"%f\n",bundle.vectorE[i]);
	}
	fclose(fp);
	
	
	///---プログラム終了処理---
	cout << endl;
	cout << "--------------" << endl;
	cout << "    Finish    " << endl;
	cout << "--------------" << endl;
	end_time_total = clock(); //実行時間計測終了
	cout << "プログラム実行時間 = " << (float)(end_time_total - start_time_total)/CLOCKS_PER_SEC << "秒" << endl << endl;


	return 0;
}



///////////////////////////////
// クイックソート用比較関数
//
int comp( const void *c1, const void *c2 )
{
	CorrPoint point1 = *(CorrPoint *)c1;
	CorrPoint point2 = *(CorrPoint *)c2;

	float tmp1 = point1.z;   /* z を基準とする */
	float tmp2 = point2.z;

	if(tmp1 == tmp2)       return 0;
	else if(tmp1 > tmp2)   return 1;
	else					 return -1;
}