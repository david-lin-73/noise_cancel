function ICAshow(sigMatrix, varargin);
% ICASHOW ---- plot the input signals with optional title and comments
%
% Command:
%    ICAshow(sigMatrix,'title','title of picture','comment','comment of the picture');
%
% Parameters:
%      'title'  ---- title of the figure
%    'comment'  ---- xlabel
%
% Author: Zhilin Zhang
%


flagcomment=0;   % ��־�ſ�ѡ����'comment'�Ƿ�ѡ�����Ϊ1��ʾ�û�ѡ��ò���������ò���ֵ 
flagtitle=0;     % ��־�ſ�ѡ����'title'�Ƿ�ѡ�����Ϊ1��ʾ�û�ѡ��ò���������ò���ֵ 
% ��ȡ��ѡ����
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs.\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'title'
                titlechar=varargin{i+1};
                flagtitle=1;
            case 'comment'
                commentchar=varargin{i+1};
                flagcomment=1;
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

rows=size(sigMatrix,1);
cols=size(sigMatrix,2);
%sigMatrix = remstd(sigMatrix);
sigMatrix = noises.ica.remstd(sigMatrix);

figure;
for i = rows : -1 : 1
    subplot(rows,1,i);
    plot(sigMatrix(i,:)); 
    
    if  flagcomment==1 & i == rows
        xlabeltext=['\bf',commentchar];
        xlabel(xlabeltext);
    end
    
    if  flagtitle==1 & i ==1
        titletext=['\bf',titlechar];
        title(titletext);
    end
     
     axis([-inf,inf,-5,5]);
end



