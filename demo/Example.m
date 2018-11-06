clear
rng('default')

N=500; % シミュレーションするステップ数
n=150; % 推定するインパルス応答のステップ数．なお，n≒NにするとFig. 1の表示が狂うことがあるが，これは最小二乗法の問題なのでdelete(hLS)とすると良い

[gtrue,u,y,sigma2,ytrue]=set_data(N,n);
% デモ用のデータセット作成関数．gtrueは真のインパルス応答，u, yは入力と観測出力，sigma2はノイズ分散，ytrueはノイズなしの出力．
% 入力に関しては白色ガウス性雑音を利用．
%% 同定部分

% 必要な情報を与えてインスタンスを生成
kernel_dataset=DCKernel(u,y,n,[]);

% ハイパーパラメータの調整．これでMyKernel内部のハイパーパラメータが変更される．
% kernel_dataset.Empirical_Bayes([1,0.8]); % 経験ベイズによるハイパーパラメータ調整．引数はハイパーパラメータの最適化の初期値．
kernel_dataset.SURE([1,0.8,0.9]); % SUREによるハイパーパラメータ調整．引数はハイパーパラメータの最適化の初期値．
% kernel_dataset.manual_hyperparameter_setting([1,0.9,0.9]); % 手動でハイパーパラメータを設定する場合はこのコマンドを利用．

[ghat,Khat]=kernel_dataset.ident(); % ghat:インパルス応答の事後平均，Khat:インパルス応答の事後分散
%% 可視化

figure(1)
clf
htrue=stairs(gtrue);
hold on
hLS=stairs(kernel_dataset.U\y,'color',[0.7,0.7,0.7]);
hhat=stairs(ghat,'r--','linewidth',1);
legend([htrue,hhat,hLS],'真のインパルス応答','カーネル法による推定値','最小二乗法による推定値')
grid on
xlim([0,n])

figure(2)
clf
hytrue=stairs(ytrue);
hold on
hy=stairs(y,'r--');
legend([hytrue,hy],'ノイズなし出力','ノイズあり出力')
