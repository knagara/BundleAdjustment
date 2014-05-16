#ifndef OPENCV_INCLUDED
#define OPENCV_INCLUDED

////プロジェクトのプロパティ⇒C\\C++⇒全般　の追加のインクルードディレクトリに
//// opencv2のあるフォルダ『C:\OpenCV\include』などを追加のこと
#include "opencv2\opencv.hpp"

using namespace cv;

#ifdef _DEBUG
    ////Debugモードの場合
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_core244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_imgproc244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_highgui244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_objdetect244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_contrib244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_features2d244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_flann244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_gpu244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_haartraining_engined.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_legacy244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_ts244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_video244d.lib")
	////追加したもの
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_ml244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_nonfree244d.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_calib3d244d.lib")

#else
    ////Releaseモードの場合
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_core244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_imgproc244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_highgui244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_objdetect244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_contrib244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_features2d244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_flann244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_gpu244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_haartraining_engined.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_legacy244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_ts244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_video244.lib")
	////追加したもの
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_ml244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_nonfree244.lib")
    #pragma comment(lib,"C:\\opencv\\build\\x86\\vc10\\lib\\opencv_calib3d244.lib")
#endif

#endif