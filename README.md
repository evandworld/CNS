# CNS
Celestial Navigation System
# 相关解释
matlab的程序，模拟飞行器（航天器）接收两个近天体和三个恒星方位时，带噪声确定自身位置。参考教材：《航天器自主天文导航原理与方法》，国防工业出版社，第二版，房建成等人著。

近天体：地球与月亮，或太阳和地球。近天体的特性是距离不是无穷大，在星敏感器中不是作为一个点，而是作为一个圆测量的，所以在测角度时，应该以其圆心为目标。
恒星：距离待测物体视为无穷大，所以只能测量角度，其在视野中的大小测量无意义。
# 建模
坐标设置：
定义航天器的坐标hhtmqi，两颗近天体的坐标：xkxk1和xkxk2.这里航天器和两个近天体用坐标表示，后面需要计算向量的时候再直接用坐标作差得到向量。

hhtmqi = [50,50,50]';
xkxk1 = [0,0,1000]'; %近天体坐标
xkxk2 = [1000,0,0]';
1
2
恒星坐标因为不受航天器位置的变化而变化（因为其足够远，所以定义航天器和近天体观测恒星时，恒星的方向矢量是不变的），所以定义三个恒星的方向矢量为s 1 s1s1、s 2 s2s2、s 3 s3s3，下面的三行代码用了矢量/矢量模的形式，得到的是单位方向矢量，个人感觉单位矢量与观测姿态角的关系像四元数与旋转矩阵的关系，但单位方向矢量理解起来似乎更直观。

s1 = [1,2,3]'/(1^2+2^2+3^2)^0.5; %恒星向量
s2 = [2,1,3]'/(2^2+1^2+3^2)^0.5;
s3 = [3,2,1]'/(3^2+2^2+1^2)^0.5;
1
2
3
与近天体1相关的计算：
计算航天器与近天体1、恒星1之间的夹角，用到两个向量的内积。航天器到近天体的矢量(定义为L11)为：航天器坐标-近天体（球心）坐标，如下式：
L 11 = h h t m q i − x k x k 1 L11 = hhtmqi - xkxk1
L11=hhtmqi−xkxk1
.
航天器与恒星1的方向矢量已经知道，是S 1 S1S1，将L 11 L11L11转化为单位矢量后，由L11、S1可以计算出来上述两个矢量之间的夹角，这个夹角即为航天器探测到的角，定义为A1，则有下式可以求得A1：
c o s ( A 1 ) = L 11 ∗ S 1 cos(A1) = L11*S1
cos(A1)=L11∗S1

A11_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s1);
A12_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s2);
A13_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s3);
1
2
3
【注】以上公式中，A11代表近地星1和恒星1之间的夹角，A12代表近地星1和恒星2之间的夹角，以此类推。A11、A12、A13合称为A1。
此处A1是标准值（理想值），在上面加上测量角度时的误差（由敏感器精度决定，最小可达1秒（1/3600度））得到实际值（测量值）：

A1 = [A11_ideal,A12_ideal,A13_ideal]'+0.001*randn(3,1);



有了夹角的测量值，则又可以通过上述公式反推得到近地星与航天器之间的方向矢量：

L1 = [cosd(A1(1)),cosd(A1(2)),cosd(A1(3))]*[s1,s2,s3]^-1; %近天体1对航天器的向量（单位向量，带噪声）


这里L1留着后面会用，接下来用同样的过程算近地星2的相关数据：

A21_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s1);
A22_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s2);
A23_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s3);
A2 = [A21_ideal,A22_ideal,A23_ideal]'+0.001*randn(3,1);
% (l1 l2 l3)(s11,s21,s31;s12,s22,s32;s13,s23,s33) = (cosa(A1),cosd(A2),cosd(A3))
L2_ideal = [cosd(A21_ideal),cosd(A22_ideal),cosd(A23_ideal)]*[s1,s2,s3]^-1;
L2 = [cosd(A2(1)),cosd(A2(2)),cosd(A2(3))]*[s1,s2,s3]^-1;

得到的L2是近地星2与航天器之间的矢量。

# 解算
由于近地星1/2的位置是知道的，所以设航天器与近地星1的距离是ρ 1 \rho 1ρ1，与近地星2的距离是ρ 2 \rho 2ρ2，则航天器与近地星1之间的矢量= L1*ρ 1 \rho1ρ1，同理，航天器与近地星2之间的矢量=L1*ρ 2 \rho2ρ2，我们是知道两个近地星的准确坐标的，所以两个矢量相减，可以得到近地星2指向近地星1的矢量，即：

L 1 ∗ ρ 1 − L 2 ∗ ρ 2 = x b x k 2 − x b x k 1 L1*\rho1-L2*\rho2=xbxk2-xbxk1
L1∗ρ1−L2∗ρ2=xbxk2−xbxk1

可列以下代码求解ρ \rhoρ:

Ac = [L1',-L2']; %vector
rho = (Ac'*Ac)^-1*Ac'*(xkxk2-xkxk1); %the length of L1 & L2(consider noisy)
1
2
代码中的ρ \rhoρ是一个矢量，里面包含了ρ 1 \rho1ρ1和ρ 2 \rho2ρ2
最后，求航天器位置，大功告成：

P1 = xkxk1+rho(1)*L1'
P2 = xkxk2+rho(2)*L2'
