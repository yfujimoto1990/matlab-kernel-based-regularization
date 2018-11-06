classdef SSKernel < kernels
    % Stable-Spline (second-order stable spline) �J�[�l���̃N���X�D
    % k(i,j) = lambda * (alpha^(i+j+max(i,j))/2-alpha^(3*max(i,j))/6)�ŗ^������J�[�l�����l����D
    % �n�C�p�[�p�����[�^theta��[lambda, alpha]�̓񎟌��x�N�g���D
    % �����lambda>0, 0<alpha<1
    
    methods
        function instance=SSKernel(u,y,n,sigma2) % �R���X�g���N�^
            % �����́C���͗�u�C�o�͗�y�C���肷��C���p���X������n�C�m�C�Y���Usigma2
            % �n�C�p�[�p�����[�^theta���������ׂ��s��������A<=b��A,b��ݒ�D
            A=[-1,0;
                   0,-1;
                   0, 1];
            b=[0;
                 0;
                 1];
            instance=instance@kernels(u,y,n,A,b,sigma2);  
        end
        
        function K=kernel_matrix(obj,theta) % ���̐ݒ�i�X�[�p�[�N���X�̃��\�b�h�j�Ƃ̊֌W��C�֐�����ύX���Ȃ����ƁD
            % �n�C�p�[�p�����[�^theta����n�~n�̃J�[�l���s��𐶐�����֐��D
            % �Ȃ��Cn���Ăяo���Ƃ���obj.n�Ƃ���Ηǂ��D
            T=ones(obj.n,1)*(1:obj.n);
            K=theta(1)*((theta(2).^(T+T'+max(T,T')))/2-(theta(2).^(3*max(T,T')))/6);
        end
       
    end
    
end