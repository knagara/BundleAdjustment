#include "stdafx.h"
#include "bundle.h"


///-----バンドル調整開始-----
void Bundle::start(selfParam _selfFlag, fixParam _fixFlag, fcParam _fcFlag, selfAlterParam _selfAlterFlag)
{
	cout << "Bundle::start" << endl;
	
	selfFlag = _selfFlag;
	selfAlterFlag = _selfAlterFlag;
	fixFlag = _fixFlag;
	fcFlag = _fcFlag;

	///パラメータ数計算
	if(selfFlag == SELF_CALIB_ON){
		cout << "SELF_CALIB_ON" << endl;
		param = 9 + paramD;
		paramK = 3 + paramD;
	}else{
		paramK = 3;
		if(fcFlag == FC_VARIABLE){
			cout << "FC_VARIABLE" << endl;
			param = 9;
		}else if(fcFlag == FC_FIX){
			cout << "FC_FIX" << endl;
			param = 6;
		}else{
			cout << "fcFlag Error" << endl;
			exit(1);
		}
	}
	allSize = 3*pointTotal + param*frameTotal;

	if(fixFlag == FIX_OFF){
		cout << "FIX_OFF" << endl;
		freedom = 0;
	}else if(fixFlag == FIX_6){
		cout << "FIX_6" << endl;
		freedom = 6;
	}else{
		cout << "FIX_7" << endl;
		freedom = 7;
	}
	useAllSize = allSize - freedom;

	///FIX7のとき，どのパラメータを固定するか
	if(fixFlag == FIX_7_AUTO){
		double minVal, maxVal;
		Point minLoc, maxLoc;
		Mat T1copy;
		T[1].copyTo(T1copy);
		T1copy.at<float>(0,0) = std::abs(T1copy.at<float>(0,0));
		T1copy.at<float>(1,0) = std::abs(T1copy.at<float>(1,0));
		T1copy.at<float>(2,0) = std::abs(T1copy.at<float>(2,0));
		minMaxLoc(T1copy, &minVal, &maxVal, &minLoc, &maxLoc);
		if(maxLoc.y == 0){
			fixFlag = FIX_7_TX;
			cout << "Fix -> tx" << endl;
		}else if(maxLoc.y == 1){
			fixFlag = FIX_7_TY;
			cout << "Fix -> ty" << endl;
		}else if(maxLoc.y == 2){
			fixFlag = FIX_7_TZ;
			cout << "Fix -> tz" << endl;
		}else{
			cout << "fixFlag Error" << endl;
			exit(1);
		}
	}

	///Mat準備
	///カメラ内部パラメータの順番は
	///tx,ty,tz,w1,w2,w3,f,cx,cy
	////一次微分 Grad : pは特徴点 fはフレームの意
	Mat Gp = (Mat_<float>(3*pointTotal,1));
	Mat Gf = (Mat_<float>(param*frameTotal,1));
	////二次微分 Hessian : Eは特徴点間 Fは特徴点とフレーム間 Gはフレーム間
	Mat HE = (Mat_<float>(3*pointTotal,3));
	Mat HF = (Mat_<float>(3*pointTotal,param*frameTotal));
	Mat HG = (Mat_<float>(param*frameTotal,param*frameTotal));
	////パラメータ固定用
	Mat HF_, HG_, Gf_;
	////パラメータ固定用 内部標定要素のみ
	Mat HG_k, Gf_k;

	///その他変数
	int bundleCount = 0; //カウンタ

	///初期再投影誤差の計算
	E = calcE();
	cout << "E = " << E << endl;
	
	vectorE.push_back(E);

	/////////////////////////////
	// 収束計算の開始
	//
	cFlag = true;
	cout << endl << "-----収束計算の開始-----" << endl;
	while(1){
		bundleCount++;
		cout << bundleCount << "回目" << endl;

		///更新値 Delta : pは特徴点 fはフレーム kは内部標定要素
		Mat Dk = (Mat_<float>(paramK*frameTotal,1));
		Dk = Scalar::all(0);

		int num;
		if(selfFlag == SELF_CALIB_ON && selfAlterFlag == ALTERNATE_ON){
			num = 6;
		}else{
			num = param;
		}

		///更新値 Delta : pは特徴点 fはフレームの意
		Mat Dp = (Mat_<float>(3*pointTotal,1));
		Mat Df = (Mat_<float>(num*frameTotal-freedom,1));
		Dp = Scalar::all(0);
		Df = Scalar::all(0);
		
		if(selfFlag == SELF_CALIB_ON && selfAlterFlag == ALTERNATE_ON){
			
			cout << "\tE    = " << E << endl;

			if(cFlag == true){
				cout << "\t内部標定要素の更新" << endl;
				///H,gの計算
				//cout << "\tcalcHg" << endl;
				calcHg(HE,HF,HG,Gp,Gf);

				HG_k = HG(Range(6*frameTotal,param*frameTotal),Range(6*frameTotal,param*frameTotal));
				Gf_k = Gf(Range(6*frameTotal,param*frameTotal),Range(0,1));
			}
			///対角成分に(1+c)をかける
			for(int i=0;i<(paramK*frameTotal);i++){
				HG_k.at<float>(i,i) *= (1+c);
			}

			///連立一次方程式の計算
			solveHg_k(HG_k,Gf_k,Dk);

			///解の更新
			update_k(Dk);

			///再投影誤差を求める
			newE = calcE();
			cout << "\tEnew = " << newE << endl;

			cout << "\t外部標定要素の更新" << endl;
			if(cFlag == true){
				///H,gの計算
				//cout << "\tcalcHg" << endl;
				calcHg(HE,HF,HG,Gp,Gf);
				
				///パラメータ固定
				if(fixFlag == FIX_OFF){
					HF_ = HF;
					HG_ = HG;
					Gf_ = Gf;
					cout << "未実装" << endl;
					exit(1);
				}else if(fixFlag == FIX_6){
					HF_ = HF(Range(0,3*pointTotal),Range(freedom,6*frameTotal));
					HG_ = HG(Range(freedom,6*frameTotal),Range(freedom,6*frameTotal));
					Gf_ = Gf(Range(freedom,6*frameTotal),Range(0,1));
				}else{
					if(fixFlag == FIX_7_TX){
						moveLine(HF,-1,-1,6,6);
						moveLine(HG,6,6,6,6);
						moveLine(Gf,6,6,-1,-1);
					}else if(fixFlag == FIX_7_TY){
						moveLine(HF,-1,-1,7,6);
						moveLine(HG,7,6,7,6);
						moveLine(Gf,7,6,-1,-1);
					}else if(fixFlag == FIX_7_TZ){
						moveLine(HF,-1,-1,8,6);
						moveLine(HG,8,6,8,6);
						moveLine(Gf,8,6,-1,-1);
					}else{
						cout << "fixFlag Error パラメータ固定時エラー" << endl;
						exit(1);
					}
					HF_ = HF(Range(0,3*pointTotal),Range(freedom,6*frameTotal));
					HG_ = HG(Range(freedom,6*frameTotal),Range(freedom,6*frameTotal));
					Gf_ = Gf(Range(freedom,6*frameTotal),Range(0,1));
				}
			}

			///対角成分に(1+c)をかける
			for(int i=0;i<pointTotal;i++){
				HE.at<float>(3*i+0,0) *= (1+c);
				HE.at<float>(3*i+1,1) *= (1+c);
				HE.at<float>(3*i+2,2) *= (1+c);
			}
			for(int i=0;i<(6*frameTotal-freedom);i++){
				HG_.at<float>(i,i) *= (1+c);
			}

			///連立一次方程式の計算
			solveHg(HE,HF_,HG_,Gp,Gf_,Dp,Df);

			///解の更新
			update_6(Dp,Df);

			///再投影誤差を求める
			newE = calcE();
			cout << "\tEnew = " << newE << endl;
				
			cout << "\te    = " << sqrt(newE/(2*KPtotal-(3*pointTotal+param*frameTotal-freedom))) << endl;

		}else{

			if(cFlag == true){
				///H,gの計算
				//cout << "\tcalcHg" << endl;
				calcHg(HE,HF,HG,Gp,Gf);

				///パラメータ固定
				if(fixFlag == FIX_OFF){
					HF_ = HF;
					HG_ = HG;
					Gf_ = Gf;
				}else if(fixFlag == FIX_6){
					HF_ = HF(Range(0,3*pointTotal),Range(freedom,param*frameTotal));
					HG_ = HG(Range(freedom,param*frameTotal),Range(freedom,param*frameTotal));
					Gf_ = Gf(Range(freedom,param*frameTotal),Range(0,1));
				}else{
					if(fixFlag == FIX_7_TX){
						moveLine(HF,-1,-1,6,6);
						moveLine(HG,6,6,6,6);
						moveLine(Gf,6,6,-1,-1);
					}else if(fixFlag == FIX_7_TY){
						moveLine(HF,-1,-1,7,6);
						moveLine(HG,7,6,7,6);
						moveLine(Gf,7,6,-1,-1);
					}else if(fixFlag == FIX_7_TZ){
						moveLine(HF,-1,-1,8,6);
						moveLine(HG,8,6,8,6);
						moveLine(Gf,8,6,-1,-1);
					}else{
						cout << "fixFlag Error パラメータ固定時エラー" << endl;
						exit(1);
					}
					HF_ = HF(Range(0,3*pointTotal),Range(freedom,param*frameTotal));
					HG_ = HG(Range(freedom,param*frameTotal),Range(freedom,param*frameTotal));
					Gf_ = Gf(Range(freedom,param*frameTotal),Range(0,1));
				}
			}

			///対角成分に(1+c)をかける
			for(int i=0;i<pointTotal;i++){
				HE.at<float>(3*i+0,0) *= (1+c);
				HE.at<float>(3*i+1,1) *= (1+c);
				HE.at<float>(3*i+2,2) *= (1+c);
			}
			for(int i=0;i<(param*frameTotal-freedom);i++){
				HG_.at<float>(i,i) *= (1+c);
			}

			///連立一次方程式の計算
			//cout << "\tsolveHg" << endl;
			clock_t start_time_solve,end_time_solve;
			start_time_solve = clock(); //実行時間計測
			solveHg(HE,HF_,HG_,Gp,Gf_,Dp,Df);
			end_time_solve = clock(); //実行時間計測終了
			cout << "計算時間 " << (float)(end_time_solve - start_time_solve)/CLOCKS_PER_SEC << "秒" << endl << endl;

			///解の更新
			//cout << "\tupdate" << endl;
			update(Dp,Df);

			///再投影誤差を求める
			newE = calcE();
			cout << "\tE    = " << E << endl;
			cout << "\tEnew = " << newE << endl;
			cout << "\te    = " << sqrt(newE/(2*KPtotal-(3*pointTotal+param*frameTotal-freedom))) << endl;
		}

		///判定
		if(newE <= E){ //Eが小さくなった
			cFlag = true;
			vectorE.push_back(newE);
			cout << "\t(E-Enew)/KPtotal = " << (E-newE)/KPtotal << endl;
			///終了条件
			///1点あたりの再投影誤差の変化量がepsilon画素以下
			if((E-newE) <= KPtotal*epsilon*epsilon){
				changeCoordToLeftUpper(); //最後に左上座標に戻す
				return; //終了
			}else{
				c *= 0.1F;
				E = newE;
			}
		}else{ //Eが大きくなってしまった
			cFlag = false;
			if(selfFlag == SELF_CALIB_ON && selfAlterFlag == ALTERNATE_ON){
				///対角成分戻す
				for(int i=0;i<(paramK*frameTotal);i++){
					HG_k.at<float>(i,i) /= (1+c);
				}
				///対角成分戻す
				for(int i=0;i<pointTotal;i++){
					HE.at<float>(3*i+0,0) /= (1+c);
					HE.at<float>(3*i+1,1) /= (1+c);
					HE.at<float>(3*i+2,2) /= (1+c);
				}
				for(int i=0;i<(6*frameTotal-freedom);i++){
					HG_.at<float>(i,i) /= (1+c);
				}
				///データ戻す
				unupdate_k(Dk);
				unupdate_6(Dp,Df);
			}else{
				///対角成分戻す
				for(int i=0;i<pointTotal;i++){
					HE.at<float>(3*i+0,0) /= (1+c);
					HE.at<float>(3*i+1,1) /= (1+c);
					HE.at<float>(3*i+2,2) /= (1+c);
				}
				for(int i=0;i<(param*frameTotal-freedom);i++){
					HG_.at<float>(i,i) /= (1+c);
				}
				///データ戻す
				unupdate(Dp,Df);
			}
			///cを10倍
			c *= 10;
			
			cout << "\tE>Enewではないので10倍します" << endl;
			cout << "\tc = " << c << endl;
			if(c > 100000000000){
				cout << "\t収束の見込みが無いため強制終了します" << endl;
				changeCoordToLeftUpper(); //最後に左上座標に戻す
				return;
			}
		}
	}

}



///-----再投影誤差の計算-----
float Bundle::calcE()
{
	//この関数だけで使用する再投影誤差E
	//return で値を返します
	float localE = 0.;

	for(int i=0;i<frameTotal;i++)
	{
		//3*4行列(I -t)を作る
		Mat T1 = (Mat_<float>(3,4) << 1, 0, 0, -T[i].at<float>(0,0), 0, 1, 0, -T[i].at<float>(1,0), 0, 0, 1, -T[i].at<float>(2,0));
		//投影行列Pを作る
		Mat Pi = K[i] * R[i].t() * T1;

		//f,cx,cy
		float f = K[i].at<float>(0,0);
		float cx = K[i].at<float>(0,2);
		float cy = K[i].at<float>(1,2);

		//Distort
		float K1 = Distort.at<float>(0,i);
		float K2 = Distort.at<float>(1,i);
		float P1 = Distort.at<float>(2,i);
		float P2 = Distort.at<float>(3,i);
		float K3 = Distort.at<float>(4,i);

		for(int j=0;j<pointTotal; j++)
		{
			if(imageX.at<float>(j,i) != (imageNaN-width/2))
			{
				//(X,Y,Z)の同次座標作る
				Mat X1 = (Mat_<float>(4,1) << points.at<float>(j,0), points.at<float>(j,1), points.at<float>(j,2), 1.);
				//p,q,r作成
				Mat Pij = Pi*X1;
				float p = Pij.at<float>(0,0);
				float q = Pij.at<float>(1,0);
				float r = Pij.at<float>(2,0);
				//x,y
				float x = p/r;
				float y = q/r;
				//x0,y0
				float x0 = imageX.at<float>(j,i);
				float y0 = imageY.at<float>(j,i);
				//Distort
				float xc = x - cx;
				float yc = y - cy;
				float r2 = xc*xc + yc*yc;
				float r4 = r2*r2;

				float r6 = r2*r2*r2;
				float dx = (K1*r2+K2*r4+K3*r6)*xc + P1*(r2+2*xc*xc) + 2*P2*xc*yc;
				float dy = (K1*r2+K2*r4+K3*r6)*yc + 2*P1*xc*yc + P2*(r2+2*yc*yc);

				///再投影誤差求めてEに加える
				//X
				localE += std::pow(p/r+dx-x0, (float)2.0);
				//Y
				localE += std::pow(q/r+dy-y0, (float)2.0);
			}
		}
	}
	return localE;
}



///-----H,gの計算-----
void Bundle::calcHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf){
	//すべてを0で初期化
	HE = Scalar::all(0);
	HF = Scalar::all(0);
	HG = Scalar::all(0);
	Gp = Scalar::all(0);
	Gf = Scalar::all(0);

	for(int i=0;i<frameTotal;i++)
	{
		//3*4行列(I -t)を作る
		Mat T1 = (Mat_<float>(3,4) << 1, 0, 0, -T[i].at<float>(0,0), 0, 1, 0, -T[i].at<float>(1,0), 0, 0, 1, -T[i].at<float>(2,0));
		//投影行列Pを作る
		Mat Pi = K[i] * R[i].t() * T1;

		//f,cx,cy
		float f = K[i].at<float>(0,0);
		float cx = K[i].at<float>(0,2);
		float cy = K[i].at<float>(1,2);

		//Distort
		float K1 = Distort.at<float>(0,i);
		float K2 = Distort.at<float>(1,i);
		float P1 = Distort.at<float>(2,i);
		float P2 = Distort.at<float>(3,i);
		float K3 = Distort.at<float>(4,i);

		///X,Y,Zに関するdp,dq,dr
		float dpX = Pi.at<float>(0,0);
		float dpY = Pi.at<float>(0,1);
		float dpZ = Pi.at<float>(0,2);
		float dqX = Pi.at<float>(1,0);
		float dqY = Pi.at<float>(1,1);
		float dqZ = Pi.at<float>(1,2);
		float drX = Pi.at<float>(2,0);
		float drY = Pi.at<float>(2,1);
		float drZ = Pi.at<float>(2,2);

		///tに関するdp,dq,dr
		Mat dpT = -(f*R[i](Range(0,3),Range(0,1))+cx*R[i](Range(0,3),Range(2,3)));
		Mat dqT = -(f*R[i](Range(0,3),Range(1,2))+cy*R[i](Range(0,3),Range(2,3)));
		Mat drT = -(R[i](Range(0,3),Range(2,3)));
		float dpTx = dpT.at<float>(0,0);
		float dpTy = dpT.at<float>(1,0);
		float dpTz = dpT.at<float>(2,0);
		float dqTx = dqT.at<float>(0,0);
		float dqTy = dqT.at<float>(1,0);
		float dqTz = dqT.at<float>(2,0);
		float drTx = drT.at<float>(0,0);
		float drTy = drT.at<float>(1,0);
		float drTz = drT.at<float>(2,0);


		for(int j=0;j<pointTotal; j++)
		{
			if(imageX.at<float>(j,i) != (imageNaN-width/2))
			{
				//(X,Y,Z)の同次座標作る
				Mat X1 = (Mat_<float>(4,1) << points.at<float>(j,0), points.at<float>(j,1), points.at<float>(j,2), 1.);
				//p,q,r
				Mat Pij = Pi*X1;
				float p = Pij.at<float>(0,0);
				float q = Pij.at<float>(1,0);
				float r = Pij.at<float>(2,0);
				//x,y
				float x = p/r;
				float y = q/r;
				//x0,y0
				float x0 = imageX.at<float>(j,i);
				float y0 = imageY.at<float>(j,i);
				//Distort
				float xc = x - cx;
				float yc = y - cy;
				float r2 = xc*xc + yc*yc;
				float r4 = r2*r2;
				float r6 = r2*r2*r2;
				float dx = (K1*r2+K2*r4+K3*r6)*xc + P1*(r2+2*xc*xc) + 2*P2*xc*yc;
				float dy = (K1*r2+K2*r4+K3*r6)*yc + 2*P1*xc*yc + P2*(r2+2*yc*yc);

				///Rに関するdp,dq,dr wを使います
				Mat X = points(Range(j,j+1),Range(0,3)).t();
				Mat dpR = (-dpT).cross(X-T[i]);
				Mat dqR = (-dqT).cross(X-T[i]);
				Mat drR = (-drT).cross(X-T[i]);
				float dpRx = dpR.at<float>(0,0);
				float dpRy = dpR.at<float>(1,0);
				float dpRz = dpR.at<float>(2,0);
				float dqRx = dqR.at<float>(0,0);
				float dqRy = dqR.at<float>(1,0);
				float dqRz = dqR.at<float>(2,0);
				float drRx = drR.at<float>(0,0);
				float drRy = drR.at<float>(1,0);
				float drRz = drR.at<float>(2,0);

				///fに関するdp,dq,dr
				float dpF = (p-cx*r)/f;
				float dqF = (q-cy*r)/f;
				float drF = 0;
				///cx,cyに関するdp,dq,dr
				float dpCx = r;
				float dqCx = 0;
				float drCx = 0;
				float dpCy = 0;
				float dqCy = r;
				float drCy = 0;
				
				///X,Y,Zに関するddx,ddy
				float ddxX = calcddx(dpX,dqX,drX,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				float ddxY = calcddx(dpY,dqY,drY,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				float ddxZ = calcddx(dpZ,dqZ,drZ,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				float ddyX = calcddy(dpX,dqX,drX,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				float ddyY = calcddy(dpY,dqY,drY,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				float ddyZ = calcddy(dpZ,dqZ,drZ,p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);

				///K1,K2,P1,P2,K3に関するddx,ddy
				float ddxK1 = r2*xc;
				float ddxK2 = r4*xc;
				float ddxK3 = r6*xc;
				float ddxP1 = r2 + 2*xc*xc;
				float ddxP2 = 2*xc*yc;
				float ddyK1 = r2*yc;
				float ddyK2 = r4*yc;
				float ddyK3 = r6*yc;
				float ddyP1 = 2*xc*yc;
				float ddyP2 = r2 + 2*yc*yc;
				
				float * dp = (float *)calloc(param, sizeof(float));
				float * dq = (float *)calloc(param, sizeof(float));
				float * dr = (float *)calloc(param, sizeof(float));
				float * ddx = (float *)calloc(param, sizeof(float));
				float * ddy = (float *)calloc(param, sizeof(float));
				int * cxFlag = (int *)calloc(param, sizeof(int));
				int * cyFlag = (int *)calloc(param, sizeof(int));

				dp[0] = dpTx;
				dp[1] = dpTy;
				dp[2] = dpTz;
				dp[3] = dpRx;
				dp[4] = dpRy;
				dp[5] = dpRz;

				dq[0] = dqTx;
				dq[1] = dqTy;
				dq[2] = dqTz;
				dq[3] = dqRx;
				dq[4] = dqRy;
				dq[5] = dqRz;

				dr[0] = drTx;
				dr[1] = drTy;
				dr[2] = drTz;
				dr[3] = drRx;
				dr[4] = drRy;
				dr[5] = drRz;

				for(int k=0;k<6;k++){
					ddx[k] = calcddx(dp[k],dq[k],dr[k],p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
					ddy[k] = calcddy(dp[k],dq[k],dr[k],p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,0,0);
				}

				///焦点距離，光軸点を更新するときのみ
				if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){

					cxFlag[7] = 1;
					cyFlag[8] = 1;

					dp[6] = dpF;
					dp[7] = dpCx;
					dp[8] = dpCy;
					dq[6] = dqF;
					dq[7] = dqCx;
					dq[8] = dqCy;
					dr[6] = drF;
					dr[7] = drCx;
					dr[8] = drCy;

					for(int k=6;k<8;k++){
						ddx[k] = calcddx(dp[k],dq[k],dr[k],p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,cxFlag[k],cyFlag[k]);
						ddy[k] = calcddy(dp[k],dq[k],dr[k],p,q,r,cx,cy,K1,K2,K3,P1,P2,xc,yc,r2,r4,r6,cxFlag[k],cyFlag[k]);
					}

					///セルフキャリブレーション
					if(selfFlag == SELF_CALIB_ON){
						ddx[9] = ddxK1;
						ddx[10] = ddxK2;
						ddx[11] = ddxP1;
						ddx[12] = ddxP2;
						ddx[13] = ddxK3;
						ddy[9] = ddyK1;
						ddy[10] = ddyK2;
						ddy[11] = ddyP1;
						ddy[12] = ddyP2;
						ddy[13] = ddyK3;
					}
				}
				///X,Y,Zの一次微分
				//dE/dXk
				Gp.at<float>(3*j+0,0) += dE(dpX,dqX,drX,p,q,r,x0,y0,dx,dy,ddxX,ddyX);
				//dE/dYk
				Gp.at<float>(3*j+1,0) += dE(dpY,dqY,drY,p,q,r,x0,y0,dx,dy,ddxY,ddyY);
				//dE/dZk
				Gp.at<float>(3*j+2,0) += dE(dpZ,dqZ,drZ,p,q,r,x0,y0,dx,dy,ddxZ,ddyZ);

				///tx,ty,txの一次微分
				//dE/dtx
				Gf.at<float>(6*i+0,0) += dE(dpTx,dqTx,drTx,p,q,r,x0,y0,dx,dy,ddx[0],ddy[0]);
				//dE/dty
				Gf.at<float>(6*i+1,0) += dE(dpTy,dqTy,drTy,p,q,r,x0,y0,dx,dy,ddx[1],ddy[1]);
				//dE/dtz
				Gf.at<float>(6*i+2,0) += dE(dpTz,dqTz,drTz,p,q,r,x0,y0,dx,dy,ddx[2],ddy[2]);

				///R: w1,w2,w3の一次微分
				//dE/dw1
				Gf.at<float>(6*i+3,0) += dE(dpRx,dqRx,drRx,p,q,r,x0,y0,dx,dy,ddx[3],ddy[3]);
				//dE/dw2
				Gf.at<float>(6*i+4,0) += dE(dpRy,dqRy,drRy,p,q,r,x0,y0,dx,dy,ddx[4],ddy[4]);
				//dE/dw3
				Gf.at<float>(6*i+5,0) += dE(dpRz,dqRz,drRz,p,q,r,x0,y0,dx,dy,ddx[5],ddy[5]);

				///焦点距離，光軸点を更新するときのみ
				if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){
					///fの一次微分
					Gf.at<float>(6*frameTotal+paramK*i+0,0) += dE(dpF,dqF,drF,p,q,r,x0,y0,dx,dy,ddx[6],ddy[6]);
					///cx,cyの一次微分
					Gf.at<float>(6*frameTotal+paramK*i+1,0) += dE(dpCx,dqCx,drCx,p,q,r,x0,y0,dx,dy,ddx[7],ddy[7]);
					Gf.at<float>(6*frameTotal+paramK*i+2,0) += dE(dpCy,dqCy,drCy,p,q,r,x0,y0,dx,dy,ddx[8],ddy[8]);
					///セルフキャリブレーション
					if(selfFlag == SELF_CALIB_ON){
						Gf.at<float>(6*frameTotal+paramK*i+3,0) += dE(0,0,0,p,q,r,x0,y0,dx,dy,ddxK1,ddyK1);
						Gf.at<float>(6*frameTotal+paramK*i+4,0) += dE(0,0,0,p,q,r,x0,y0,dx,dy,ddxK2,ddyK2);
						Gf.at<float>(6*frameTotal+paramK*i+5,0) += dE(0,0,0,p,q,r,x0,y0,dx,dy,ddxP1,ddyP1);
						Gf.at<float>(6*frameTotal+paramK*i+6,0) += dE(0,0,0,p,q,r,x0,y0,dx,dy,ddxP2,ddyP2);
						Gf.at<float>(6*frameTotal+paramK*i+7,0) += dE(0,0,0,p,q,r,x0,y0,dx,dy,ddxK3,ddyK3);
					}
				}


				///X,Y,Z同士の二次微分
				///上三角半分のみ あとで下三角半分にコピーする

				//dXdX
				HE.at<float>(3*j+0,0) += ddE(dpX,dqX,drX,dpX,dqX,drX,p,q,r,ddxX,ddyX,ddxX,ddyX);
				//dXdY
				HE.at<float>(3*j+0,1) += ddE(dpX,dqX,drX,dpY,dqY,drY,p,q,r,ddxX,ddyX,ddxY,ddyY);
				//dXdZ
				HE.at<float>(3*j+0,2) += ddE(dpX,dqX,drX,dpZ,dqZ,drZ,p,q,r,ddxX,ddyX,ddxZ,ddyZ);
				//dYdY
				HE.at<float>(3*j+1,1) += ddE(dpY,dqY,drY,dpY,dqY,drY,p,q,r,ddxY,ddyY,ddxY,ddyY);
				//dYdZ
				HE.at<float>(3*j+1,2) += ddE(dpY,dqY,drY,dpZ,dqZ,drZ,p,q,r,ddxY,ddyY,ddxZ,ddyZ);
				//dZdZ
				HE.at<float>(3*j+2,2) += ddE(dpZ,dqZ,drZ,dpZ,dqZ,drZ,p,q,r,ddxZ,ddyZ,ddxZ,ddyZ);


				///フレーム同士の二次微分
				///上三角半分のみ あとで下三角半分にコピーする

				///カメラ同士：並進・回転の二次微分
				for(int k=0;k<6;k++){
					for(int l=0;l<k+1;l++){
						HG.at<float>(6*i+l,6*i+k) += ddE(dp[l],dq[l],dr[l],dp[k],dq[k],dr[k],p,q,r,ddx[l],ddy[l],ddx[k],ddy[k]);
					}
				}
				///特徴点とカメラ：並進・回転の二次微分
				for(int k=0;k<6;k++){
					HF.at<float>(3*j+0,6*i+k) += ddE(dpX,dqX,drX,dp[k],dq[k],dr[k],p,q,r,ddxX,ddyX,ddx[k],ddy[k]);
					HF.at<float>(3*j+1,6*i+k) += ddE(dpY,dqY,drY,dp[k],dq[k],dr[k],p,q,r,ddxY,ddyY,ddx[k],ddy[k]);
					HF.at<float>(3*j+2,6*i+k) += ddE(dpZ,dqZ,drZ,dp[k],dq[k],dr[k],p,q,r,ddxZ,ddyZ,ddx[k],ddy[k]);
				}

				
				///焦点距離，光軸点を更新するときのみ
				if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){
					///カメラ同士：焦点距離・光軸点の二次微分
					for(int k=0;k<paramK;k++){
						for(int l=0;l<k+1;l++){
							HG.at<float>(6*frameTotal+paramK*i+l,6*frameTotal+paramK*i+k) += ddE(dp[l+6],dq[l+6],dr[l+6],dp[k+6],dq[k+6],dr[k+6],p,q,r,ddx[l+6],ddy[l+6],ddx[k+6],ddy[k+6]);
						}
					}
					///カメラ同士：焦点距離・光軸点と並進・回転の二次微分
					for(int k=0;k<paramK;k++){
						for(int l=0;l<6;l++){
							HG.at<float>(6*i+l,6*frameTotal+paramK*i+k) += ddE(dp[l],dq[l],dr[l],dp[k+6],dq[k+6],dr[k+6],p,q,r,ddx[l],ddy[l],ddx[k+6],ddy[k+6]);
						}
					}
					///特徴点とカメラ：焦点距離・光軸点の二次微分
					for(int k=0;k<paramK;k++){
						HF.at<float>(3*j+0,6*frameTotal+paramK*i+k) += ddE(dpX,dqX,drX,dp[k+6],dq[k+6],dr[k+6],p,q,r,ddxX,ddyX,ddx[k+6],ddy[k+6]);
						HF.at<float>(3*j+1,6*frameTotal+paramK*i+k) += ddE(dpY,dqY,drY,dp[k+6],dq[k+6],dr[k+6],p,q,r,ddxY,ddyY,ddx[k+6],ddy[k+6]);
						HF.at<float>(3*j+2,6*frameTotal+paramK*i+k) += ddE(dpZ,dqZ,drZ,dp[k+6],dq[k+6],dr[k+6],p,q,r,ddxZ,ddyZ,ddx[k+6],ddy[k+6]);
					}
				}
			}
		}
	}

	///二次微分の下三角半分コピー
	for(int j=0;j<pointTotal; j++){
		HE.at<float>(3*j+1,0) = HE.at<float>(3*j+0,1);
		HE.at<float>(3*j+2,0) = HE.at<float>(3*j+0,2);
		HE.at<float>(3*j+2,1) = HE.at<float>(3*j+1,2);
	}
	for(int i=0;i<frameTotal;i++){
		for(int j=0;j<6;j++){
			for(int k=0;k<j;k++){
				HG.at<float>(6*i+j,6*i+k) = HG.at<float>(6*i+k,6*i+j);
			}
		}
		///焦点距離，光軸点を更新するときのみ
		if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){
			for(int j=0;j<paramK;j++){
				for(int k=0;k<j;k++){
					HG.at<float>(6*frameTotal+paramK*i+j,6*frameTotal+paramK*i+k) = HG.at<float>(6*frameTotal+paramK*i+k,6*frameTotal+paramK*i+j);
				}
			}
			for(int j=0;j<paramK;j++){
				for(int k=0;k<6;k++){
					HG.at<float>(6*frameTotal+paramK*i+j,6*i+k) = HG.at<float>(6*i+k,6*frameTotal+paramK*i+j);
				}
			}
		}
	}
}




///-----連立一次方程式の計算-----
void Bundle::solveHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf,Mat& Dp,Mat& Df)
{
	//cout << "\t  Dfの計算" << endl;

	int num;
	if(selfFlag == SELF_CALIB_ON && selfAlterFlag == ALTERNATE_ON){
		num = 6;
	}else{
		num = param;
	}
	
	///左の係数部分 G-FEF
	Mat FEF = (Mat_<float>(num*frameTotal-freedom,num*frameTotal-freedom));
	FEF = Scalar::all(0);
	///右の係数部分 FEGp-Gf
	Mat FEGp = (Mat_<float>(num*frameTotal-freedom,1));
	FEGp = Scalar::all(0);

	for(int i=0;i<pointTotal;i++){
		Mat Einv = HE(Range(3*i,3*i+3),Range(0,3)).inv();
		Mat Fi = HF(Range(3*i,3*i+3),Range(0,num*frameTotal-freedom));
		Mat Gpi = Gp(Range(3*i,3*i+3),Range(0,1));

		FEF +=  Fi.t() * Einv * Fi;
		FEGp += Fi.t() * Einv * Gpi;
	}

	///連立一次方程式を解いてDfを求める
	cv::solve(HG-FEF, FEGp-Gf, Df);


	//cout << "\t  Dpの計算" << endl;
	for(int i=0;i<pointTotal;i++){
		//Mat用意
		Mat Dpi = (Mat_<float>(3,1));
		Mat Einv = HE(Range(3*i,3*i+3),Range(0,3)).inv();
		Mat Fi = HF(Range(3*i,3*i+3),Range(0,num*frameTotal-freedom));
		Mat Gpi = Gp(Range(3*i,3*i+3),Range(0,1));
		//Dpi計算
		Dpi = -(Einv*(Fi*Df+Gpi));
		//代入
		Dp.at<float>(3*i+0,0) = Dpi.at<float>(0,0);
		Dp.at<float>(3*i+1,0) = Dpi.at<float>(1,0);
		Dp.at<float>(3*i+2,0) = Dpi.at<float>(2,0);
	}
}





///-----連立一次方程式の計算-----
void Bundle::solveHg_k(Mat& HG,Mat& Gf,Mat& Dk)
{
	for(int i=0;i<frameTotal;i++){
		Mat HGinv = HG(Range(paramK*i,paramK*(i+1)),Range(paramK*i,paramK*(i+1))).inv();
		Mat Gfi = Gf(Range(paramK*i,paramK*(i+1)),Range(0,1));
		Mat Dki = HGinv * Gfi;
		//代入
		Dk.at<float>(paramK*i+0,0) = Dki.at<float>(0,0);
		Dk.at<float>(paramK*i+1,0) = Dki.at<float>(1,0);
		Dk.at<float>(paramK*i+2,0) = Dki.at<float>(2,0);
		Dk.at<float>(paramK*i+3,0) = Dki.at<float>(3,0);
		Dk.at<float>(paramK*i+4,0) = Dki.at<float>(4,0);
		Dk.at<float>(paramK*i+5,0) = Dki.at<float>(5,0);
		Dk.at<float>(paramK*i+6,0) = Dki.at<float>(6,0);
		Dk.at<float>(paramK*i+7,0) = Dki.at<float>(7,0);
	}
}




///-----解の更新-----
void Bundle::update_k(Mat& Dk)
{
	for(int i=0;i<frameTotal;i++){
		K[i].at<float>(0,0) += Dk.at<float>(paramK*i+0,0);
		K[i].at<float>(1,1) += Dk.at<float>(paramK*i+0,0);
		K[i].at<float>(0,2) += Dk.at<float>(paramK*i+1,0);
		K[i].at<float>(1,2) += Dk.at<float>(paramK*i+2,0);
		Distort.at<float>(0,i) += Dk.at<float>(paramK*i+3,0);
		Distort.at<float>(1,i) += Dk.at<float>(paramK*i+4,0);
		Distort.at<float>(2,i) += Dk.at<float>(paramK*i+5,0);
		Distort.at<float>(3,i) += Dk.at<float>(paramK*i+6,0);
		Distort.at<float>(4,i) += Dk.at<float>(paramK*i+7,0);
	}
}

///-----解を戻す-----
void Bundle::unupdate_k(Mat& Dk)
{
	for(int i=0;i<frameTotal;i++){
		K[i].at<float>(0,0) -= Dk.at<float>(paramK*i+0,0);
		K[i].at<float>(1,1) -= Dk.at<float>(paramK*i+0,0);
		K[i].at<float>(0,2) -= Dk.at<float>(paramK*i+1,0);
		K[i].at<float>(1,2) -= Dk.at<float>(paramK*i+2,0);
		Distort.at<float>(0,i) -= Dk.at<float>(paramK*i+3,0);
		Distort.at<float>(1,i) -= Dk.at<float>(paramK*i+4,0);
		Distort.at<float>(2,i) -= Dk.at<float>(paramK*i+5,0);
		Distort.at<float>(3,i) -= Dk.at<float>(paramK*i+6,0);
		Distort.at<float>(4,i) -= Dk.at<float>(paramK*i+7,0);
	}
}





///-----解の更新-----
void Bundle::update_6(Mat& Dp,Mat& Df)
{
	///3次元位置の更新
	for(int i=0;i<pointTotal;i++){
		points.at<float>(i,0) += Dp.at<float>(3*i+0,0);
		points.at<float>(i,1) += Dp.at<float>(3*i+1,0);
		points.at<float>(i,2) += Dp.at<float>(3*i+2,0);
	}

	///カメラ変数の更新
	if(fixFlag == FIX_7_TX){
		T[1].at<float>(1,0) += Df.at<float>(0,0);
		T[1].at<float>(2,0) += Df.at<float>(1,0);
	}else if(fixFlag == FIX_7_TY){
		T[1].at<float>(0,0) += Df.at<float>(0,0);
		T[1].at<float>(2,0) += Df.at<float>(1,0);
	}else if(fixFlag == FIX_7_TZ){
		T[1].at<float>(0,0) += Df.at<float>(0,0);
		T[1].at<float>(1,0) += Df.at<float>(1,0);
	}
	Mat Rodri = (Mat_<float>(3,3));
	Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
	R[1] = Rodri*R[1];

	for(int i=0;i<frameTotal-2;i++){
		//並進
		T[i+2].at<float>(0,0) += Df.at<float>(5+6*i+0,0);
		T[i+2].at<float>(1,0) += Df.at<float>(5+6*i+1,0);
		T[i+2].at<float>(2,0) += Df.at<float>(5+6*i+2,0);
		//回転
		Mat Rodri = (Mat_<float>(3,3));
		Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
		R[i+2] = Rodri*R[i+2];
	}
}





///-----解を戻す-----
void Bundle::unupdate_6(Mat& Dp,Mat& Df)
{
	///3次元位置の更新
	for(int i=0;i<pointTotal;i++){
		points.at<float>(i,0) -= Dp.at<float>(3*i+0,0);
		points.at<float>(i,1) -= Dp.at<float>(3*i+1,0);
		points.at<float>(i,2) -= Dp.at<float>(3*i+2,0);
	}

	///カメラ変数の更新
	if(fixFlag == FIX_7_TX){
		T[1].at<float>(1,0) -= Df.at<float>(0,0);
		T[1].at<float>(2,0) -= Df.at<float>(1,0);
	}else if(fixFlag == FIX_7_TY){
		T[1].at<float>(0,0) -= Df.at<float>(0,0);
		T[1].at<float>(2,0) -= Df.at<float>(1,0);
	}else if(fixFlag == FIX_7_TZ){
		T[1].at<float>(0,0) -= Df.at<float>(0,0);
		T[1].at<float>(1,0) -= Df.at<float>(1,0);
	}
	Mat Rodri = (Mat_<float>(3,3));
	Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
	R[1] = Rodri.t()*R[1];

	for(int i=0;i<frameTotal-2;i++){
		//並進
		T[i+2].at<float>(0,0) -= Df.at<float>(5+6*i+0,0);
		T[i+2].at<float>(1,0) -= Df.at<float>(5+6*i+1,0);
		T[i+2].at<float>(2,0) -= Df.at<float>(5+6*i+2,0);
		//回転
		Mat Rodri = (Mat_<float>(3,3));
		Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
		R[i+2] = Rodri.t()*R[i+2];
	}
}




///-----解の更新-----
void Bundle::update(Mat& Dp,Mat& Df)
{
	///3次元位置の更新
	for(int i=0;i<pointTotal;i++){
		points.at<float>(i,0) += Dp.at<float>(3*i+0,0);
		points.at<float>(i,1) += Dp.at<float>(3*i+1,0);
		points.at<float>(i,2) += Dp.at<float>(3*i+2,0);
	}

	///カメラ変数の更新
	////FIX_OFF////
	if(fixFlag == FIX_OFF){
		if(selfFlag == SELF_CALIB_ON){
			cout << "実装してません" << endl;
			exit(1);
		}
		for(int i=0;i<frameTotal;i++){
			//並進
			T[i].at<float>(0,0) += Df.at<float>(6*i+0,0);
			T[i].at<float>(1,0) += Df.at<float>(6*i+1,0);
			T[i].at<float>(2,0) += Df.at<float>(6*i+2,0);
			//回転
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
			R[i] = Rodri*R[i];
			///焦点距離，光軸点を更新するときのみ
			if(fcFlag == FC_VARIABLE){
				K[i].at<float>(0,0) += Df.at<float>(6*frameTotal+3*i+0,0);
				K[i].at<float>(1,1) += Df.at<float>(6*frameTotal+3*i+0,0);
				K[i].at<float>(0,2) += Df.at<float>(6*frameTotal+3*i+1,0);
				K[i].at<float>(1,2) += Df.at<float>(6*frameTotal+3*i+2,0);
			}
		}
	////FIX_6////
	}else if(fixFlag == FIX_6){
		if(selfFlag == SELF_CALIB_ON){
			cout << "実装してません" << endl;
			exit(1);
		}
		///焦点距離，光軸点を更新するときのみ
		if(fcFlag == FC_VARIABLE){
			K[0].at<float>(0,0) += Df.at<float>(6*(frameTotal-1)+0,0);
			K[0].at<float>(1,1) += Df.at<float>(6*(frameTotal-1)+0,0);
			K[0].at<float>(0,2) += Df.at<float>(6*(frameTotal-1)+1,0);
			K[0].at<float>(1,2) += Df.at<float>(6*(frameTotal-1)+2,0);
			for(int i=0;i<frameTotal-1;i++){
				//並進
				T[i+1].at<float>(0,0) += Df.at<float>(6*i+0,0);
				T[i+1].at<float>(1,0) += Df.at<float>(6*i+1,0);
				T[i+1].at<float>(2,0) += Df.at<float>(6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
				R[i+1] = Rodri*R[i+1];
				K[i+1].at<float>(0,0) += Df.at<float>(6*(frameTotal-1)+3*i+0,0);
				K[i+1].at<float>(1,1) += Df.at<float>(6*(frameTotal-1)+3*i+0,0);
				K[i+1].at<float>(0,2) += Df.at<float>(6*(frameTotal-1)+3*i+1,0);
				K[i+1].at<float>(1,2) += Df.at<float>(6*(frameTotal-1)+3*i+2,0);
			}
		}
		///焦点距離，光軸点を更新しない
		else if(fcFlag == FC_FIX){
			for(int i=0;i<frameTotal-1;i++){
				//並進
				T[i+1].at<float>(0,0) += Df.at<float>(6*i+0,0);
				T[i+1].at<float>(1,0) += Df.at<float>(6*i+1,0);
				T[i+1].at<float>(2,0) += Df.at<float>(6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
				R[i+1] = Rodri*R[i+1];
			}
		}else{
			cout << "fcFlag Error at update" << endl;
			exit(1);
		}
	////FIX_7////
	}else{
		///焦点距離，光軸点を更新するときのみ
		if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){
			K[0].at<float>(0,0) += Df.at<float>(6*frameTotal-7+0,0);
			K[0].at<float>(1,1) += Df.at<float>(6*frameTotal-7+0,0);
			K[0].at<float>(0,2) += Df.at<float>(6*frameTotal-7+1,0);
			K[0].at<float>(1,2) += Df.at<float>(6*frameTotal-7+2,0);
			if(selfFlag == SELF_CALIB_ON){
				Distort.at<float>(0,0) += Df.at<float>(6*frameTotal-7+3,0);
				Distort.at<float>(1,0) += Df.at<float>(6*frameTotal-7+4,0);
				Distort.at<float>(2,0) += Df.at<float>(6*frameTotal-7+5,0);
				Distort.at<float>(3,0) += Df.at<float>(6*frameTotal-7+6,0);
				Distort.at<float>(4,0) += Df.at<float>(6*frameTotal-7+7,0);
			}
			if(fixFlag == FIX_7_TX){
				T[1].at<float>(1,0) += Df.at<float>(0,0);
				T[1].at<float>(2,0) += Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TY){
				T[1].at<float>(0,0) += Df.at<float>(0,0);
				T[1].at<float>(2,0) += Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TZ){
				T[1].at<float>(0,0) += Df.at<float>(0,0);
				T[1].at<float>(1,0) += Df.at<float>(1,0);
			}
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
			R[1] = Rodri*R[1];
			K[1].at<float>(0,0) += Df.at<float>(6*frameTotal-7+paramK+0,0);
			K[1].at<float>(1,1) += Df.at<float>(6*frameTotal-7+paramK+0,0);
			K[1].at<float>(0,2) += Df.at<float>(6*frameTotal-7+paramK+1,0);
			K[1].at<float>(1,2) += Df.at<float>(6*frameTotal-7+paramK+2,0);
			if(selfFlag == SELF_CALIB_ON){
				Distort.at<float>(0,1) += Df.at<float>(6*frameTotal-7+paramK+3,0);
				Distort.at<float>(1,1) += Df.at<float>(6*frameTotal-7+paramK+4,0);
				Distort.at<float>(2,1) += Df.at<float>(6*frameTotal-7+paramK+5,0);
				Distort.at<float>(3,1) += Df.at<float>(6*frameTotal-7+paramK+6,0);
				Distort.at<float>(4,1) += Df.at<float>(6*frameTotal-7+paramK+7,0);
			}

			for(int i=0;i<frameTotal-2;i++){
				//並進
				T[i+2].at<float>(0,0) += Df.at<float>(5+6*i+0,0);
				T[i+2].at<float>(1,0) += Df.at<float>(5+6*i+1,0);
				T[i+2].at<float>(2,0) += Df.at<float>(5+6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
				R[i+2] = Rodri*R[i+2];
				K[i+2].at<float>(0,0) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+0,0);
				K[i+2].at<float>(1,1) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+0,0);
				K[i+2].at<float>(0,2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+1,0);
				K[i+2].at<float>(1,2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+2,0);
				if(selfFlag == SELF_CALIB_ON){
					Distort.at<float>(0,i+2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+3,0);
					Distort.at<float>(1,i+2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+4,0);
					Distort.at<float>(2,i+2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+5,0);
					Distort.at<float>(3,i+2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+6,0);
					Distort.at<float>(4,i+2) += Df.at<float>(6*frameTotal-7+paramK*(i+2)+7,0);
				}
			}
		}
		///焦点距離，光軸点を更新しない
		else if(fcFlag == FC_FIX){	
			if(fixFlag == FIX_7_TX){
				T[1].at<float>(1,0) += Df.at<float>(0,0);
				T[1].at<float>(2,0) += Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TY){
				T[1].at<float>(0,0) += Df.at<float>(0,0);
				T[1].at<float>(2,0) += Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TZ){
				T[1].at<float>(0,0) += Df.at<float>(0,0);
				T[1].at<float>(1,0) += Df.at<float>(1,0);
			}
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
			R[1] = Rodri*R[1];

			for(int i=0;i<frameTotal-2;i++){
				//並進
				T[i+2].at<float>(0,0) += Df.at<float>(5+6*i+0,0);
				T[i+2].at<float>(1,0) += Df.at<float>(5+6*i+1,0);
				T[i+2].at<float>(2,0) += Df.at<float>(5+6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
				R[i+2] = Rodri*R[i+2];
			}
		}else{
			cout << "fcFlag Error at update" << endl;
			exit(1);
		}
	}
}




///-----解を戻す-----
void Bundle::unupdate(Mat& Dp,Mat& Df)
{
	///3次元位置の更新
	for(int i=0;i<pointTotal;i++){
		points.at<float>(i,0) -= Dp.at<float>(3*i+0,0);
		points.at<float>(i,1) -= Dp.at<float>(3*i+1,0);
		points.at<float>(i,2) -= Dp.at<float>(3*i+2,0);
	}

	///カメラ変数の更新
	////FIX_OFF////
	if(fixFlag == FIX_OFF){
		if(selfFlag == SELF_CALIB_ON){
			cout << "実装してません" << endl;
			exit(1);
		}
		for(int i=0;i<frameTotal;i++){
			//並進
			T[i].at<float>(0,0) -= Df.at<float>(6*i+0,0);
			T[i].at<float>(1,0) -= Df.at<float>(6*i+1,0);
			T[i].at<float>(2,0) -= Df.at<float>(6*i+2,0);
			//回転
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
			R[i] = Rodri.t()*R[i];
			///焦点距離，光軸点を更新するときのみ
			if(fcFlag == FC_VARIABLE){
				K[i].at<float>(0,0) -= Df.at<float>(6*frameTotal+3*i+0,0);
				K[i].at<float>(1,1) -= Df.at<float>(6*frameTotal+3*i+0,0);
				K[i].at<float>(0,2) -= Df.at<float>(6*frameTotal+3*i+1,0);
				K[i].at<float>(1,2) -= Df.at<float>(6*frameTotal+3*i+2,0);
			}
		}
	////FIX_6////
	}else if(fixFlag == FIX_6){
		if(selfFlag == SELF_CALIB_ON){
			cout << "実装してません" << endl;
			exit(1);
		}
		///焦点距離，光軸点を更新するときのみ
		if(fcFlag == FC_VARIABLE){
			K[0].at<float>(0,0) -= Df.at<float>(6*(frameTotal-1)+0,0);
			K[0].at<float>(1,1) -= Df.at<float>(6*(frameTotal-1)+0,0);
			K[0].at<float>(0,2) -= Df.at<float>(6*(frameTotal-1)+1,0);
			K[0].at<float>(1,2) -= Df.at<float>(6*(frameTotal-1)+2,0);
			for(int i=0;i<frameTotal-1;i++){
				//並進
				T[i+1].at<float>(0,0) -= Df.at<float>(6*i+0,0);
				T[i+1].at<float>(1,0) -= Df.at<float>(6*i+1,0);
				T[i+1].at<float>(2,0) -= Df.at<float>(6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
				R[i+1] = Rodri.t()*R[i+1];
				K[i+1].at<float>(0,0) -= Df.at<float>(6*(frameTotal-1)+3*i+0,0);
				K[i+1].at<float>(1,1) -= Df.at<float>(6*(frameTotal-1)+3*i+0,0);
				K[i+1].at<float>(0,2) -= Df.at<float>(6*(frameTotal-1)+3*i+1,0);
				K[i+1].at<float>(1,2) -= Df.at<float>(6*(frameTotal-1)+3*i+2,0);
			}
		}
		///焦点距離，光軸点を更新しない
		else if(fcFlag == FC_FIX){
			for(int i=0;i<frameTotal-1;i++){
				//並進
				T[i+1].at<float>(0,0) -= Df.at<float>(6*i+0,0);
				T[i+1].at<float>(1,0) -= Df.at<float>(6*i+1,0);
				T[i+1].at<float>(2,0) -= Df.at<float>(6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(6*i+3,6*i+6),Range(0,1)),Rodri);
				R[i+1] = Rodri.t()*R[i+1];
			}
		}else{
			cout << "fcFlag Error at update" << endl;
			exit(1);
		}
	////FIX_7////
	}else{
		///焦点距離，光軸点を更新するときのみ
		if(fcFlag == FC_VARIABLE || selfFlag == SELF_CALIB_ON){
			K[0].at<float>(0,0) -= Df.at<float>(6*frameTotal-7+0,0);
			K[0].at<float>(1,1) -= Df.at<float>(6*frameTotal-7+0,0);
			K[0].at<float>(0,2) -= Df.at<float>(6*frameTotal-7+1,0);
			K[0].at<float>(1,2) -= Df.at<float>(6*frameTotal-7+2,0);
			if(selfFlag == SELF_CALIB_ON){
				Distort.at<float>(0,0) -= Df.at<float>(6*frameTotal-7+3,0);
				Distort.at<float>(1,0) -= Df.at<float>(6*frameTotal-7+4,0);
				Distort.at<float>(2,0) -= Df.at<float>(6*frameTotal-7+5,0);
				Distort.at<float>(3,0) -= Df.at<float>(6*frameTotal-7+6,0);
				Distort.at<float>(4,0) -= Df.at<float>(6*frameTotal-7+7,0);
			}
			if(fixFlag == FIX_7_TX){
				T[1].at<float>(1,0) -= Df.at<float>(0,0);
				T[1].at<float>(2,0) -= Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TY){
				T[1].at<float>(0,0) -= Df.at<float>(0,0);
				T[1].at<float>(2,0) -= Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TZ){
				T[1].at<float>(0,0) -= Df.at<float>(0,0);
				T[1].at<float>(1,0) -= Df.at<float>(1,0);
			}
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
			R[1] = Rodri.t()*R[1];
			K[1].at<float>(0,0) -= Df.at<float>(6*frameTotal-7+paramK+0,0);
			K[1].at<float>(1,1) -= Df.at<float>(6*frameTotal-7+paramK+0,0);
			K[1].at<float>(0,2) -= Df.at<float>(6*frameTotal-7+paramK+1,0);
			K[1].at<float>(1,2) -= Df.at<float>(6*frameTotal-7+paramK+2,0);
			if(selfFlag == SELF_CALIB_ON){
				Distort.at<float>(0,1) -= Df.at<float>(6*frameTotal-7+paramK+3,0);
				Distort.at<float>(1,1) -= Df.at<float>(6*frameTotal-7+paramK+4,0);
				Distort.at<float>(2,1) -= Df.at<float>(6*frameTotal-7+paramK+5,0);
				Distort.at<float>(3,1) -= Df.at<float>(6*frameTotal-7+paramK+6,0);
				Distort.at<float>(4,1) -= Df.at<float>(6*frameTotal-7+paramK+7,0);
			}

			for(int i=0;i<frameTotal-2;i++){
				//並進
				T[i+2].at<float>(0,0) -= Df.at<float>(5+6*i+0,0);
				T[i+2].at<float>(1,0) -= Df.at<float>(5+6*i+1,0);
				T[i+2].at<float>(2,0) -= Df.at<float>(5+6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
				R[i+2] = Rodri.t()*R[i+2];
				K[i+2].at<float>(0,0) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+0,0);
				K[i+2].at<float>(1,1) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+0,0);
				K[i+2].at<float>(0,2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+1,0);
				K[i+2].at<float>(1,2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+2,0);
				if(selfFlag == SELF_CALIB_ON){
					Distort.at<float>(0,i+2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+3,0);
					Distort.at<float>(1,i+2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+4,0);
					Distort.at<float>(2,i+2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+5,0);
					Distort.at<float>(3,i+2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+6,0);
					Distort.at<float>(4,i+2) -= Df.at<float>(6*frameTotal-7+paramK*(i+2)+7,0);
				}
			}
		}
		///焦点距離，光軸点を更新しない
		else if(fcFlag == FC_FIX){	
			if(fixFlag == FIX_7_TX){
				T[1].at<float>(1,0) -= Df.at<float>(0,0);
				T[1].at<float>(2,0) -= Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TY){
				T[1].at<float>(0,0) -= Df.at<float>(0,0);
				T[1].at<float>(2,0) -= Df.at<float>(1,0);
			}else if(fixFlag == FIX_7_TZ){
				T[1].at<float>(0,0) -= Df.at<float>(0,0);
				T[1].at<float>(1,0) -= Df.at<float>(1,0);
			}
			Mat Rodri = (Mat_<float>(3,3));
			Rodrigues(Df(Range(2,5),Range(0,1)),Rodri);
			R[1] = Rodri.t()*R[1];

			for(int i=0;i<frameTotal-2;i++){
				//並進
				T[i+2].at<float>(0,0) -= Df.at<float>(5+6*i+0,0);
				T[i+2].at<float>(1,0) -= Df.at<float>(5+6*i+1,0);
				T[i+2].at<float>(2,0) -= Df.at<float>(5+6*i+2,0);
				//回転
				Mat Rodri = (Mat_<float>(3,3));
				Rodrigues(Df(Range(5+6*i+3,5+6*i+6),Range(0,1)),Rodri);
				R[i+2] = Rodri.t()*R[i+2];
			}
		}else{
			cout << "fcFlag Error at update" << endl;
			exit(1);
		}
	}
}



///-----特定の行と列を移動させる-----
void Bundle::moveLine(Mat& A,int _row,int row,int _col,int col)
{
	///row
	if(_row == -1 || row == -1 || _row == row){
		//何もしない
	}else{
		Mat Row;
		A.row(_row).copyTo(Row);
		if(row > _row){
			for(int i=_row+1;i<=row;i++){
				for(int j=0;j<A.cols;j++){
					A.at<float>(i-1,j) = A.at<float>(i,j);
				}
			}
		}else{
			for(int i=_row-1;i>=row;i--){
				for(int j=0;j<A.cols;j++){
					A.at<float>(i+1,j) = A.at<float>(i,j);
				}
			}
		}
		for(int j=0;j<A.cols;j++){
			A.at<float>(row,j) = Row.at<float>(0,j);
		}
	}

	///col
	if(_col == -1 || col == -1 || _col == col){
		//何もしない
	}else{
		Mat Col;
		A.col(_col).copyTo(Col);
		if(col > _col){
			for(int i=_col+1;i<=col;i++){
				for(int j=0;j<A.rows;j++){
					A.at<float>(j,i-1) = A.at<float>(j,i);
				}
			}
		}else{
			for(int i=_col-1;i>=col;i--){
				for(int j=0;j<A.rows;j++){
					A.at<float>(j,i+1) = A.at<float>(j,i);
				}
			}
		}
		for(int j=0;j<A.rows;j++){
			A.at<float>(j,col) = Col.at<float>(j,0);
		}
	}
}



///-----左上座標から中心座標に変換-----
void Bundle::changeCoordToCenter()
{
	//画像座標
	for(int i=0;i<pointTotal;i++){
		for(int j=0;j<frameTotal;j++){
			imageX.at<float>(i,j) -= width/2;
			imageY.at<float>(i,j) -= height/2;
		}
	}
	//カメラ主点位置
	for(int i=0;i<frameTotal;i++){
		K[i].at<float>(0,2) -= width/2;
		K[i].at<float>(1,2) -= height/2;
	}
}


///-----中心座標から左上座標に変換-----
void Bundle::changeCoordToLeftUpper()
{
	//画像座標
	for(int i=0;i<pointTotal;i++){
		for(int j=0;j<frameTotal;j++){
			imageX.at<float>(i,j) += width/2;
			imageY.at<float>(i,j) += height/2;
		}
	}
	//カメラ主点位置
	for(int i=0;i<frameTotal;i++){
		K[i].at<float>(0,2) += width/2;
		K[i].at<float>(1,2) += height/2;
	}
}


///-----焦点距離をfx,fyからfに統一-----
void Bundle::fxfyAverage()
{
	for(int i=0;i<frameTotal;i++){
		float fx = K[i].at<float>(0,0);
		float fy = K[i].at<float>(1,1);
		K[i].at<float>(0,0) = (fx+fy)/2;
		K[i].at<float>(1,1) = (fx+fy)/2;
	}
}

///-----KPtotalの計算-----
void Bundle::calcKPtotal()
{
	KPtotal = 0;
	for(int i=0;i<pointTotal;i++){
		for(int j=0;j<frameTotal;j++){
			if(imageX.at<float>(i,j) != imageNaN - width/2){
				KPtotal++;
			}
		}
	}
}