classdef DCKernel < kernels
    % Diagonal-Correlated  カーネルのクラス．
    % k(i,j) = lambda * alpha^((i+j)/2)*rho^(abs(i-j))で与えられるカーネルを考える．
    % ハイパーパラメータthetaは[lambda, alpha,rho]の三次元ベクトル．
    % 制約はlambda>0, 0<alpha<1, -1<rho<1
    
    methods
        function instance=DCKernel(u,y,n,sigma2) % コンストラクタ
            % 引数は，入力列u，出力列y，推定するインパルス応答長n，ノイズ分散sigma2
            % ハイパーパラメータthetaが満たすべき不等式制約A<=bのA,bを設定．
            A=[-1, 0, 0;
                   0,-1, 0;
                   0, 1, 0;
                   0, 0,-1;
                   0, 0, 1];
            b=[0;
                 0;
                 1;
                 1;
                 1];
            instance=instance@kernels(u,y,n,A,b,sigma2);  
        end
        
        function K=kernel_matrix(obj,theta) % 他の設定（スーパークラスのメソッド）との関係上，関数名を変更しないこと．
            % ハイパーパラメータthetaからn×nのカーネル行列を生成する関数．
            % なお，nを呼び出すときはobj.nとすれば良い．
            T=ones(obj.n,1)*(1:obj.n);
            K=theta(1)*(theta(2).^((T+T')/2)).*(theta(3).^(abs(T-T')));
        end
       
    end
    
end