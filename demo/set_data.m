function [gtrue,u,y,sigma2,ytrue]=set_data(N,n)

% �ΏۃV�X�e���ݒ�Fz�̈�ł̑��Ύ�����1�ȏ�ɂȂ�悤�ɂ��Ă��������D
% ���Ύ�����1���傫���Ƃނ����Ԃ�����CTC�J�[�l���Ȃǂł͏��������̐��萸�x�������邱�Ƃ�����܂��D
sys=zpk(0,[0.7,0.8],1,[]);

gtrue=impulse(sys,1:n); % �^�̃C���p���X����

% �m�C�Y���U
sigma2=10;

%% �ϑ��f�[�^���W
u=randn(N,1); % ���́F���F���K�K�E�X�G��

% �o�͌v�Z�Flsim�ł͒��B�������Ԃ����C��X�̘g�g�݂ł͒��B���͂Ȃ����Ƃ͊��m�Ƃ��Ă���̂Œ��B�����͏ȗ�
ytrue=lsim(sys,[u;0]);
ytrue=ytrue(2:end);

% �m�C�Y�̕t��
y=ytrue+sqrt(sigma2)*randn(N,1);

end