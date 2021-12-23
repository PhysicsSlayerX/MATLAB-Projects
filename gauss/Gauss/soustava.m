function soustava;
global lambda R1 W1 scrsz q1 a1 b1 k;

handle=findobj('Tag','vlndelka');
lambda=str2num(get(handle,'String'));
handle=findobj('Tag','sirka');
R1=str2num(get(handle,'String'));
handle=findobj('Tag','polomer');
W1=str2num(get(handle,'String'));
handle=findobj('Tag','Vyber');
W1=W1*1e-3;
vyber=get(handle,'Value');
close;
lambda=lambda*1e-9;
k=2*pi/lambda;
a1=k^2*R1*W1^4/(k^2*W1^4+4*R1^2);
b1=-2*k*R1^2*W1^2/(k^2*W1^2+4*R1^2);
q1=a1+j*b1;
if vyber == 1
   prostor;
elseif vyber == 2
   cocka;
elseif vyber == 3
   rov_zrcadlo;   
elseif vyber == 4
   kul_zrcadlo;
elseif vyber == 5
   rov_rozhrani;   
elseif vyber == 6
   kul_rozhrani;   
end
