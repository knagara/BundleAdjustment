BundleAdjustment
================
  
####Self-calibration Bundle Adjustment  
  
[INPUT]  
* Data set of corresponding points（(X,Y,Z) in global coordinate，(x,y) on image coordinate on each frame）  
* Translation and orientation (R,T) on each frame  
* Interior parameters (focal length = f, principal point = (cx, cy))  on each frame 
  
[OUTPUT]  
* Point cloud（(X,Y,Z) in global coordinate）  
* Translation and orientation (R,T) on each frame  
* Interior parameters (focal length = f, principal point = (cx, cy))  on each frame  
* Interior parameters (lens distortion parameters = (K1,K2,K3,P1,P2))  on each frame   
  
================
  
####セルフキャリブレーション付きバンドル調整  
  
入力  
　対応点データ（x,y,z空間座標，x,y画像座標）  
　カメラ外部標定要素（並進，回転）  
　カメラ内部標定要素（焦点距離，主点位置）  
  
出力  
　点群座標（x,y,z空間座標）  
　カメラ外部標定要素（並進，回転）  
　カメラ内部標定要素（焦点距離，主点位置）  
　カメラ内部標定要素（レンズ歪み係数）  
  
  
