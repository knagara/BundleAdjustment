#ifndef BUNDLE_H_INCLUDED
#define BUNDLE_H_INCLUDED

#include "stdafx.h"
#include "opencv.h"


//セルフキャリブレーションするかどうか
enum selfParam{
	SELF_CALIB_ON = 0,
	SELF_CALIB_OFF = 1
};

//交互セルフキャリブレーション
enum selfAlterParam{
	ALTERNATE_ON = 0,
	ALTERNATE_OFF = 1
};

//どのようにパラメータを固定するか
enum fixParam{
	//7つ固定，R0,T0とあともう一つ
	FIX_7_AUTO = 0, //自動判定
	FIX_7_TX = 1, //txを固定
	FIX_7_TY = 2, //txを固定
	FIX_7_TZ = 3, //txを固定
	//6つ固定，R0,T0
	FIX_6 = 4,
	//固定しない。パラメータを全て動かす
	FIX_OFF = 5
};

//カメラの焦点距離，主点位置を更新するかどうか
enum fcParam{
	FC_VARIABLE = 0,
	FC_FIX = 1
};

class Bundle
{
private:
	///private変数
	Mat& points; //点群(X,Y,Z)
	Mat& imageX; //画像座標x
	Mat& imageY; //画像座標y
	Mat * R; //回転行列 3*3
	Mat * T; //並進ベクトル 3*1
	Mat * K; //カメラ内部行列 3*3

	//点の総数，フレーム数，画像サイズ
	const int pointTotal, frameTotal, width, height;

	//空間上の点がカメラに投影されてないときの画像座標
	//コンストラクタで-999.に初期化しています
	float imageNaN;

	//レーベンバーグ・マーカート法のc
	//コンストラクタで0.00001Fに初期化しています
	float c;

	//終了条件のイプシロン
	//コンストラクタで0.1に初期化しています
	float epsilon; 

	//歪み係数の数
	//コンストラクタで5に初期化しています
	//4 : K1,K2,P1,P2
	//5 : K1,K2,K3,P1,P2
	int paramD;

public:
	///public変数
	selfParam selfFlag; //セルフキャリブレーションするかどうか
	selfAlterParam selfAlterFlag; //交互セルフキャリブレーション
	fixParam fixFlag; //どのようにパラメータを固定するか
	fcParam fcFlag; //カメラの焦点距離，主点位置を更新するかどうか
	float E, newE; //再投影誤差E
	int allSize; //全パラメータ
	int useAllSize; //更新するパラメータ
	int freedom; //自由度
	int param; //カメラのパラメータ数 6or10
	int paramK; //カメラの内部パラメータ数 3 or (3+paramD)
	bool cFlag; //反復計算のcのフラグ
	int KPtotal; //各フレームに写っている特徴点の数を合計したもの
	Mat Distort; //歪み係数
	vector<float> vectorE; //再投影誤差Eを保存するvector


	///-----コンストラクタ-----
	//private変数に，参照渡しされた引数が代入されます
	//呼び出し元のデータをそのまま使うので，このクラス内で値が変更されると呼び出し元のデータの値も変更されます。
	Bundle(Mat& _points, Mat& _imageX, Mat& _imageY, Mat * _R, Mat * _T, Mat * _K, int _width, int _height) : points(_points), imageX(_imageX), imageY(_imageY), width(_width), height(_height), imageNaN(-999.), epsilon(0.1F), c(0.0001F), pointTotal(imageX.rows), frameTotal(imageX.cols), paramD(5)
	{
		R = _R;
		T = _T;
		K = _K;

		//歪み係数の初期化
		Distort = (Mat_<float>(paramD,frameTotal));
		Distort = Scalar::all(0);

		//左上座標から中心座標に変換
		changeCoordToCenter();

		//焦点距離をfx,fyからfに統一
		fxfyAverage();

		//KPtotalの計算
		calcKPtotal();
		
		cout << "Bundle::Constructor" << endl;
		cout << "frameTotal = " << frameTotal << endl;
		cout << "pointTotal = " << pointTotal << endl;
		cout << "KPTotal = " << KPtotal << endl << endl;

	}

	///-----バンドル調整開始-----
	void start(selfParam _selfFlag=SELF_CALIB_ON, fixParam _fixFlag=FIX_7_AUTO, fcParam _fcFlag=FC_VARIABLE, selfAlterParam _selfAlterFlag=ALTERNATE_ON);

	///-----再投影誤差Eの計算-----
	float calcE();

	///-----H,gの計算-----
	void calcHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf);

	///-----連立一次方程式の計算-----
	void solveHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf,Mat& Dp,Mat& Df);

	///-----連立一次方程式の計算-----
	void solveHg_k(Mat& HG,Mat& Gf,Mat& Dk);

	///-----解の更新-----
	void update(Mat& Dp,Mat& Df);

	///-----解の更新-----
	void update_6(Mat& Dp,Mat& Df);

	///-----解の更新-----
	void update_k(Mat& Dk);

	///-----解を戻す-----
	void unupdate(Mat& Dp,Mat& Df);

	///-----解を戻す-----
	void unupdate_6(Mat& Dp,Mat& Df);

	///-----解を戻す-----
	void unupdate_k(Mat& Dk);

	///-----特定の行と列を移動させる-----
	void moveLine(Mat& A,int i_=-1,int i=-1,int j_=-1,int j=-1);

	///-----左上座標から中心座標に変換-----
	void changeCoordToCenter();

	///-----中心座標から左上座標に変換-----
	void changeCoordToLeftUpper();

	///-----焦点距離をfx,fyからfに統一-----
	void fxfyAverage();

	///-----KPtotalの計算-----
	void calcKPtotal();

	///-----imageNaN変更-----
	//imageNaNの初期値は-999.です
	void set_imageNaN(float _imageNaN)
	{
		cout << "Set imageNaN " << imageNaN;
		imageNaN = _imageNaN;
		cout << " to " << imageNaN << endl << endl;
	}

	///-----c変更-----
	//cの初期値は0.00001Fです
	void set_c(float _c)
	{
		cout << "Set c " << c;
		c = _c;
		cout << " to " << c << endl << endl;
	}

	///-----epsilon変更-----
	//epsilonの初期値は0.1です
	void set_epsilon(float _epsilon)
	{
		cout << "Set epsilon " << epsilon;
		epsilon = _epsilon;
		cout << " to " << epsilon << endl << endl;
	}

	///-----paramD変更-----
	//paramDの初期値は5です
	void set_paramD(int _paramD)
	{
		cout << "Set paramD " << paramD;
		paramD = _paramD;
		cout << " to " << paramD << endl << endl;

		//歪み係数の初期化
		Distort = (Mat_<float>(paramD,frameTotal));
		Distort = Scalar::all(0);
	}

	///-----Distortデータ入力-----
	void set_Distort(Mat& _Distort)
	{
		cout << "Set Distort " << endl;

		Distort = _Distort;
	}

	///-----歪み係数の微分-----
	inline float calcddx(float dp, float dq, float dr, float p, float q, float r, float cx, float cy, float K1, float K2, float K3, float P1, float P2, float xc, float yc, float r2, float r4, float r6, int cxFlag, int cyFlag)
	{
		float x_ = (r*dp-p*dr)/(r*r);
		float y_ = (r*dq-q*dr)/(r*r);
		float r2_; 
		if(cxFlag==0 && cyFlag==0){
			r2_ = 2*xc*x_ + 2*yc*y_;
		}else if(cxFlag!=0){
			r2_ = 2*xc*(x_-1) + 2*yc*y_;
		}else if(cyFlag!=0){
			r2_ = 2*xc*x_ + 2*yc*(y_-1);
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}
		float r4_ = 2*r2*r2_;
		float r6_ = 3*r4*r2_;

		float ddx;
		if(cxFlag==0 && cyFlag==0){
			ddx = (K1*r2+K2*r4+K3*r6)*x_ + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*x_) + 2*P2*(x_*yc+xc*y_);
		}else if(cxFlag!=0){
			ddx = (K1*r2+K2*r4+K3*r6)*(x_-1) + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*(x_-1)) + 2*P2*((x_-1)*yc+xc*y_);
		}else if(cyFlag!=0){
			ddx = (K1*r2+K2*r4+K3*r6)*x_ + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*x_) + 2*P2*(x_*yc+xc*(y_-1));
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}

		return ddx;
	}

	///-----歪み係数の微分-----
	inline float calcddy(float dp, float dq, float dr, float p, float q, float r, float cx, float cy, float K1, float K2, float K3, float P1, float P2, float xc, float yc, float r2, float r4, float r6, int cxFlag, int cyFlag)
	{
		float x_ = (r*dp-p*dr)/(r*r);
		float y_ = (r*dq-q*dr)/(r*r);
		float r2_; 
		if(cxFlag==0 && cyFlag==0){
			r2_ = 2*xc*x_ + 2*yc*y_;
		}else if(cxFlag!=0){
			r2_ = 2*xc*(x_-1) + 2*yc*y_;
		}else if(cyFlag!=0){
			r2_ = 2*xc*x_ + 2*yc*(y_-1);
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}
		float r4_ = 2*r2*r2_;
		float r6_ = 3*r4*r2_;

		float ddy;
		if(cxFlag==0 && cyFlag==0){
			ddy = (K1*r2+K2*r4+K3*r6)*y_ + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*(x_*yc+xc*y_) + P2*(r2_+4*yc*y_);
		}else if(cxFlag!=0){
			ddy = (K1*r2+K2*r4+K3*r6)*y_ + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*((x_-1)*yc+xc*y_) + P2*(r2_+4*yc*y_);
		}else if(cyFlag!=0){
			ddy = (K1*r2+K2*r4+K3*r6)*(y_-1) + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*(x_*yc+xc*(y_-1)) + P2*(r2_+4*yc*(y_-1));
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}

		return ddy;
	}
	
	///-----一次微分の計算-----
	inline float dE(float dp, float dq, float dr, float p, float q, float r, float x0, float y0, float dx=0, float dy=0, float ddx=0, float ddy=0)
	{
		return 2*((p/r+dx-x0)*((r*dp-p*dr)/(r*r)+ddx)+(q/r+dy-y0)*((r*dq-q*dr)/(r*r)+ddy));
	}
	
	///-----二次微分の計算-----
	inline float ddE(float dpk, float dqk, float drk, float dpl, float dql, float drl, float p, float q, float r, float ddxk=0, float ddyk=0, float ddxl=0, float ddyl=0)
	{
		return 2*(((r*dpk-p*drk)/(r*r)+ddxk)*((r*dpl-p*drl)/(r*r)+ddxl)+((r*dqk-q*drk)/(r*r)+ddyk)*((r*dql-q*drl)/(r*r)+ddyl));
	}

};

#endif