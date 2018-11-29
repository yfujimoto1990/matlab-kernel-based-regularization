classdef SSKernel < kernels
    % Stable-Spline (second-order stable spline) カーネルのクラス．
    % k(i,j) = lambda * (alpha^(i+j+max(i,j))/2-alpha^(3*max(i,j))/6)で与えられるカーネルを考える．
    % ハイパーパラメータthetaは[lambda, alpha]の二次元ベクトル．
    % 制約はlambda>0, 0<alpha<1
    
    methods
        function instance=SSKernel(u,y,n,sigma2) % コンストラクタ
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
            K=theta(1)*((theta(2).^(T+T'+max(T,T')))/2-(theta(2).^(3*max(T,T')))/6);
        end
       
    end
    
end
