# matlab-kernel-based-regularization
カーネルを用いたシステム同定のMATLABでの実装サンプルです．
2018/11/29 matlab 2016b以前にも対応できるよう，使用していたスーパークラスを変更しました．


使用する際にはoptimization toolboxが必要です．

新たに設計したカーネルを容易に実装できるよう，オブジェクト指向でプログラミングしました．

Example.m：デモ用のプログラムサンプル
kernels.m：カーネルを用いたシステム同定のためのクラス．他のクラスファイルはこれを継承している
TCkernel.m：Tuned-Correlated kernel (1st order stable spline kernel)のクラス
DCkernel.m：Diagonal-Correlated kernelのクラス
SSkernel.m：(2nd order) Stable-Spline kernelのクラス
set_data.m：デモで用いるデータを生成する関数


参考文献：
G. Pillonetto et al. "Kernel methods in system identification, machine learning and function estimation: A survey." Automatica, Vol. 50, No. 3, pp. 657-682 (2014).
B. Mu, T. Chen, and L. Ljung. "On asymptotic properties of hyperparameter estimators for kernel-based regularization methods." Automatica, Vol. 94, pp.381-395 (2018).


