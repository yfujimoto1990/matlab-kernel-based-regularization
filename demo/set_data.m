function [gtrue,u,y,sigma2,ytrue]=set_data(N,n)

% 対象システム設定：z領域での相対次数が1以上になるようにしてください．
% 相対次数が1より大きいとむだ時間があり，TCカーネルなどでは初期応答の推定精度が落ちることがあります．
sys=zpk(0,[0.7,0.8],1,[]);

gtrue=impulse(sys,1:n); % 真のインパルス応答

% ノイズ分散
sigma2=10;

%% 観測データ収集
u=randn(N,1); % 入力：白色正規ガウス雑音

% 出力計算：lsimでは直達部分も返すが，我々の枠組みでは直達項はないことは既知としているので直達部分は省略
ytrue=lsim(sys,[u;0]);
ytrue=ytrue(2:end);

% ノイズの付加
y=ytrue+sqrt(sigma2)*randn(N,1);

end