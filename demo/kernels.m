classdef kernels < matlab.mixin.SetGet 
    % カーネル全般に関わるパラメータやメソッドを規定するクラス．
    % ハイパーパラメータの決定，与えられたデータからの同定などはこちら．
    
    properties 
        theta % ハイパーパラメータ
        constraint_A % ハイパーパラメータの制約：A*theta<=bなる制約のAに対応
        constraint_b % ハイパーパラメータの制約：A*theta<=bなる制約のbに対応         
        U % 入力列の成すToeplitz行列．インパルス応答gに対し，出力がUgで与えられるように構成する
        y % 観測出力ベクトル
        Nfir % 推定するインパルス応答の長さ
        sigma2 % ノイズ分散
        N % 観測データ数
        K % カーネル行列．n×nのサイズ
    end
    
    
    methods    
    
        function instance=kernels(u,y,n,A,b,sigma2) % コンストラクタ
            instance.Nfir=n;
            instance.y=reshape(y,[],1); % yとして横ベクトルが入ってきたとしても縦に変換
            instance.N=length(y);
            instance.U=instance.toeplitz_for_convolution(u,n);
            if ~isempty(sigma2)
                instance.sigma2=sigma2;
            else
                gLS=instance.U\y;
                instance.sigma2=var(y-instance.U*(gLS));
            end
            instance.constraint_A=A;
            instance.constraint_b=b;
        end
                
        function [ghat,Khat]=ident(obj)
            % 与えられたデータと設定されたハイパーパラメータを用いてインパルス応答を同定
            if isempty(obj.theta)
                % ハイパーパラメータが設定されていない場合エラーを返す
                error('先にハイパーパラメータを決定してください')
            end
            
            obj.K=obj.kernel_matrix(obj.theta); % 事前共分散を計算
            
            % 事後共分散の計算．事前共分散のランクに応じて計算を変更．
            if rank(obj.K)<obj.Nfir 
                % ランク落ちしている場合の計算.逆行列を用いないがKhatが（数値的に）半正定対称行列にならない場合が存在する．
                Khat=obj.K-obj.K*obj.U'*((obj.U*obj.K*obj.U'+obj.sigma2*eye(obj.N))\obj.U*obj.K);
            else % 事前共分散がフルランクなら，数値的にはこちらの方が安定している．
                Khat=inv(inv(obj.K)+obj.U'*obj.U/obj.sigma2);
            end
            
            ghat=obj.K*((obj.sigma2*eye(obj.Nfir)+obj.U'*obj.U*obj.K)\(obj.U'*obj.y)); % 事後平均を計算．
        end
        
        function manual_hyperparameter_setting(obj,theta)
            % 手動でハイパーパラメータを設定するメソッド．
            obj.theta=reshape(theta,[],1);            
        end
        
        function Empirical_Bayes(obj,x0)
            % 周辺尤度最大化によるハイパーパラメータ調整．xがハイパーパラメータに対応．
            Z0=@(x) obj.U*obj.kernel_matrix(x)*obj.U'+obj.sigma2*eye(obj.N); % 出力の事前共分散
            Z=@(x) (Z0(x)+Z0(x)')/2; % 数値的に共分散行列が対称行列とならない場合があるので，その補正
            
            J=@(x) sum(log(eig(Z(x))))+obj.y'*(Z(x)\obj.y); % 対数周辺尤度関数．定数項は除いている
            
            obj.theta=fmincon(J,x0,obj.constraint_A,obj.constraint_b); % fminconを用いた最適化．
        end
        
        function SURE(obj,x0)
            % SUREによるハイパーパラメータ調整．xがハイパーパラメータに対応
            ghat_theta=@(x)  obj.kernel_matrix(x)*((obj.sigma2*eye(obj.Nfir)+(obj.U'*obj.U)*obj.kernel_matrix(x))\(obj.U'*obj.y));% パラメータxの元での推定インパルス応答
            Jsure=@(x) norm(obj.y-obj.U*ghat_theta(x))^2+2*obj.sigma2*trace(obj.U*obj.kernel_matrix(x)*((obj.U'*obj.U*obj.kernel_matrix(x)+obj.sigma2*eye(obj.Nfir))\obj.U'));
            
            obj.theta=fmincon(Jsure,x0,obj.constraint_A,obj.constraint_b); % fminconを用いた最適化．
        end
        
        function value=get.constraint_A(obj)
            value=obj.constraint_A;
        end
        function value=get.constraint_b(obj)
            value=obj.constraint_b;
        end
        
        function set.theta(obj,theta)
            theta=reshape(theta,[],1); %　横ベクトルで入ってきても対応
            % thetaの次元数が正しく設定されているか判定．
            % その後不等式制約を満たしているか確認
            A=obj.get.constraint_A();
            b=obj.get.constraint_b(); 
            ntheta=size(A,2);
            if length(theta)==ntheta
                if nnz(A*theta<=b)==length(b)
                    obj.theta=theta;
                else
                    error('ハイパーパラメータが制約を満たしていません')
                end
            else
                error('ハイパーパラメータは%d 次元のベクトルとしてください\n', ntheta)
            end
        end
        
    end

    methods(Static) 
        function U=toeplitz_for_convolution(u,n) 
            % 入力ベクトルUとインパルス応答長nから適切なUを構成する．Uについてはプロパティのコメント参照
            U=toeplitz(u,[u(1),zeros(1,n-1)]);            
        end
    end    
    
end
