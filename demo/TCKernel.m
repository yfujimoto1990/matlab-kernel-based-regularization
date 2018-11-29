classdef TCKernel < kernels
    % Tuned-Correlated (first-order stable spline) カーネルのクラス．
    % k(i,j) = lambda * alpha^(max(i,j))で与えられるカーネルを考える．
    % ハイパーパラメータthetaは[lambda, alpha]の二次元ベクトル．
    % 制約はlambda>0, 0<alpha<1
    
    methods
        function instance=TCKernel(u,y,n,sigma2) % コンストラクタ
            % 引数は，入力列u，出力列y，推定するインパルス応答長n，ノイズ分散sigma2
            % ハイパーパラメータthetaが満たすべき不等式制約A<=bのA,bを設定．
            A=[-1,0;
                   0,-1;
                   0, 1];
            b=[0;
                 0;
                 1];
            instance=instance@kernels(u,y,n,A,b,sigma2);  
        end
        
        function K=kernel_matrix(obj,theta) % 他の設定（スーパークラスのメソッド）との関係上，関数名を変更しないこと．
            % ハイパーパラメータthetaからn×nのカーネル行列を生成する関数．
            % なお，nを呼び出すときはobj.Nfirとすれば良い．
            T=ones(obj.Nfir,1)*(1:obj.Nfir);
            K=theta(1)*theta(2).^max(T,T');
        end
       
    end
    
end
