clear
rng('default')

N=500; % �V�~�����[�V��������X�e�b�v��
n=150; % ���肷��C���p���X�����̃X�e�b�v���D�Ȃ��Cn��N�ɂ����Fig. 1�̕\�����������Ƃ����邪�C����͍ŏ����@�̖��Ȃ̂�delete(hLS)�Ƃ���Ɨǂ�

[gtrue,u,y,sigma2,ytrue]=set_data(N,n);
% �f���p�̃f�[�^�Z�b�g�쐬�֐��Dgtrue�͐^�̃C���p���X�����Cu, y�͓��͂Ɗϑ��o�́Csigma2�̓m�C�Y���U�Cytrue�̓m�C�Y�Ȃ��̏o�́D
% ���͂Ɋւ��Ă͔��F�K�E�X���G���𗘗p�D
%% ���蕔��

% �K�v�ȏ���^���ăC���X�^���X�𐶐�
kernel_dataset=DCKernel(u,y,n,[]);

% �n�C�p�[�p�����[�^�̒����D�����MyKernel�����̃n�C�p�[�p�����[�^���ύX�����D
% kernel_dataset.Empirical_Bayes([1,0.8]); % �o���x�C�Y�ɂ��n�C�p�[�p�����[�^�����D�����̓n�C�p�[�p�����[�^�̍œK���̏����l�D
kernel_dataset.SURE([1,0.8,0.9]); % SURE�ɂ��n�C�p�[�p�����[�^�����D�����̓n�C�p�[�p�����[�^�̍œK���̏����l�D
% kernel_dataset.manual_hyperparameter_setting([1,0.9,0.9]); % �蓮�Ńn�C�p�[�p�����[�^��ݒ肷��ꍇ�͂��̃R�}���h�𗘗p�D

[ghat,Khat]=kernel_dataset.ident(); % ghat:�C���p���X�����̎��㕽�ρCKhat:�C���p���X�����̎��㕪�U
%% ����

figure(1)
clf
htrue=stairs(gtrue);
hold on
hLS=stairs(kernel_dataset.U\y,'color',[0.7,0.7,0.7]);
hhat=stairs(ghat,'r--','linewidth',1);
legend([htrue,hhat,hLS],'�^�̃C���p���X����','�J�[�l���@�ɂ�鐄��l','�ŏ����@�ɂ�鐄��l')
grid on
xlim([0,n])

figure(2)
clf
hytrue=stairs(ytrue);
hold on
hy=stairs(y,'r--');
legend([hytrue,hy],'�m�C�Y�Ȃ��o��','�m�C�Y����o��')
