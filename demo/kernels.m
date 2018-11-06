classdef kernels < matlab.mixin.SetGetExactNames 
    % �J�[�l���S�ʂɊւ��p�����[�^�⃁�\�b�h���K�肷��N���X�D
    % �n�C�p�[�p�����[�^�̌���C�^����ꂽ�f�[�^����̓���Ȃǂ͂�����D
    
    properties 
        theta % �n�C�p�[�p�����[�^
        constraint_A % �n�C�p�[�p�����[�^�̐���FA*theta<=b�Ȃ鐧���A�ɑΉ�
        constraint_b % �n�C�p�[�p�����[�^�̐���FA*theta<=b�Ȃ鐧���b�ɑΉ�         
        U % ���͗�̐���Toeplitz�s��D�C���p���X����g�ɑ΂��C�o�͂�Ug�ŗ^������悤�ɍ\������
        y % �ϑ��o�̓x�N�g��
        n % ���肷��C���p���X�����̒���
        sigma2 % �m�C�Y���U
        N % �ϑ��f�[�^��
        K % �J�[�l���s��Dn�~n�̃T�C�Y
    end
    
    
    methods    
    
        function instance=kernels(u,y,n,A,b,sigma2) % �R���X�g���N�^
            instance.n=n;
            instance.y=reshape(y,[],1); % y�Ƃ��ĉ��x�N�g���������Ă����Ƃ��Ă��c�ɕϊ�
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
            % �^����ꂽ�f�[�^�Ɛݒ肳�ꂽ�n�C�p�[�p�����[�^��p���ăC���p���X�����𓯒�
            if isempty(obj.theta)
                % �n�C�p�[�p�����[�^���ݒ肳��Ă��Ȃ��ꍇ�G���[��Ԃ�
                error('��Ƀn�C�p�[�p�����[�^�����肵�Ă�������')
            end
            
            obj.K=obj.kernel_matrix(obj.theta); % ���O�����U���v�Z
            
            % ���㋤���U�̌v�Z�D���O�����U�̃����N�ɉ����Čv�Z��ύX�D
            if rank(obj.K)<obj.n 
                % �����N�������Ă���ꍇ�̌v�Z.�t�s���p���Ȃ���Khat���i���l�I�Ɂj������Ώ̍s��ɂȂ�Ȃ��ꍇ�����݂���D
                Khat=obj.K-obj.K*obj.U'*((obj.U*obj.K*obj.U'+obj.sigma2*eye(obj.N))\obj.U*obj.K);
            else % ���O�����U���t�������N�Ȃ�C���l�I�ɂ͂�����̕������肵�Ă���D
                Khat=inv(inv(obj.K)+obj.U'*obj.U/obj.sigma2);
            end
            
            ghat=obj.K*((obj.sigma2*eye(obj.n)+obj.U'*obj.U*obj.K)\(obj.U'*obj.y)); % ���㕽�ς��v�Z�D
        end
        
        function manual_hyperparameter_setting(obj,theta)
            % �蓮�Ńn�C�p�[�p�����[�^��ݒ肷�郁�\�b�h�D
            obj.theta=reshape(theta,[],1);            
        end
        
        function Empirical_Bayes(obj,x0)
            % ���Ӗޓx�ő剻�ɂ��n�C�p�[�p�����[�^�����Dx���n�C�p�[�p�����[�^�ɑΉ��D
            Z0=@(x) obj.U*obj.kernel_matrix(x)*obj.U'+obj.sigma2*eye(obj.N); % �o�͂̎��O�����U
            Z=@(x) (Z0(x)+Z0(x)')/2; % ���l�I�ɋ����U�s�񂪑Ώ̍s��ƂȂ�Ȃ��ꍇ������̂ŁC���̕␳
            
            J=@(x) sum(log(eig(Z(x))))+obj.y'*(Z(x)\obj.y); % �ΐ����Ӗޓx�֐��D�萔���͏����Ă���
            
            obj.theta=fmincon(J,x0,obj.constraint_A,obj.constraint_b); % fmincon��p�����œK���D
        end
        
        function SURE(obj,x0)
            % SURE�ɂ��n�C�p�[�p�����[�^�����Dx���n�C�p�[�p�����[�^�ɑΉ�
            ghat_theta=@(x)  obj.kernel_matrix(x)*((obj.sigma2*eye(obj.n)+(obj.U'*obj.U)*obj.kernel_matrix(x))\(obj.U'*obj.y));% �p�����[�^x�̌��ł̐���C���p���X����
            Jsure=@(x) norm(obj.y-obj.U*ghat_theta(x))^2+2*obj.sigma2*trace(obj.U*obj.kernel_matrix(x)*((obj.U'*obj.U*obj.kernel_matrix(x)+obj.sigma2*eye(obj.n))\obj.U'));
            
            obj.theta=fmincon(Jsure,x0,obj.constraint_A,obj.constraint_b); % fmincon��p�����œK���D
        end
        
        function value=get.constraint_A(obj)
            value=obj.constraint_A;
        end
        function value=get.constraint_b(obj)
            value=obj.constraint_b;
        end
        
        function set.theta(obj,theta)
            theta=reshape(theta,[],1); %�@���x�N�g���œ����Ă��Ă��Ή�
            % theta�̎��������������ݒ肳��Ă��邩����D
            % ���̌�s��������𖞂����Ă��邩�m�F
            A=obj.get.constraint_A();
            b=obj.get.constraint_b(); 
            ntheta=size(A,2);
            if length(theta)==ntheta
                if nnz(A*theta<=b)==length(b)
                    obj.theta=theta;
                else
                    error('�n�C�p�[�p�����[�^������𖞂����Ă��܂���')
                end
            else
                error('�n�C�p�[�p�����[�^��%d �����̃x�N�g���Ƃ��Ă�������\n', ntheta)
            end
        end
        
    end

    methods(Static) 
        function U=toeplitz_for_convolution(u,n) 
            % ���̓x�N�g��U�ƃC���p���X������n����K�؂�U���\������DU�ɂ��Ă̓v���p�e�B�̃R�����g�Q��
            U=toeplitz(u,[u(1),zeros(1,n-1)]);            
        end
    end    
    
end