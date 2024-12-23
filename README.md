# mmsMD_R0
--- 
1、命名规则  
2、源文件功能说明  
3、输入文件规则说明
--- 
## 命名规则
仿原版程序，t_开头代表type，h_开头代表host变量等。纯cpu程序也使用host变量  
原程序用的para关键词有点烦人，parameter和parallel都是para，难以区分，本程序要尽可能避免  
本程序中，para仅用于parameter，而parallel不缩写  
nb仅用于neighbor，nonbonded缩写成nonb  
## 源文件功能说明
### main.cpp
仅起到一个启动器的作用，不想让这个文件太复杂，主体放在了md.cpp上
argv用于接收终端传来的信息，例如
./BRMD pretreat做选择的预处理
./BRMD 不加其他词条，就是按默认方式跑MD
### includes.h
给头文件打个包，免得每个文件开头一大串
### md.cpp
t_nstep 管理各种输入输出操作的频率（几步一做）  
参数初始化init_all
单步循环onestep
build_lj_para，预先计算一部分12-6势的参数
### dsa.h
1.模板类无序数组array；特点是动态扩增，可以删除中间一个元素并将末尾元素摘下补上。对于不需要排序查找，需要插入删除而又频繁读写的数据（如结构信息），可保证全部操作O(n)复杂度以下。
array.reduce():归约求和函数，在速度要求低精度要求高的位置最好采用这个，可以避免大数吃小数
### file.cpp
1.ini读写系列函数，给出配置文件的键值，返回数据；如key(键值)=1000(数据)  
其核心是strtok，根据=号将字符串切割；再用atoi等将字符串转换为数字  
现版本，等号左右可随意加空格

2.

### top.h
核心数据，原子的速度和位置。这里有两种处理方式，即结构体里装array<速度>，array<位置>和直接构造array<结构体（速度、位置）>两种方案。这里要结合内存读写速度和安全性选择。原程序t_atomprop选择了前者，为了方便后续移植我们直接继承这一思路。原程序用了atoms描述质量电荷等静态数据、atomprop描述运动学数据，不如我们改个名字吧，更清晰：t_atom_kinematics和t_atom_statics  
为了方便，把nstep(从ini获得的输出操作步数)和radius(各种截断半径)
### nblist.cpp
这里用了一个struct来构造nblist，毕竟后续是要改造的。  
二维数组的方式也是权宜之计，这个结构在GPU并行计算的时候可能会遇到麻烦
周期性边界条件的处理抄自原程序：  
dx -= bx * rintf(dx * bxinv);  
rintf来自math.h(C99新特性)将浮点数的小数部分四舍五入得到一个浮点数.
### bond.cpp/angle.cpp
思路相同。读取molecule列表，从molecule对象里获得分子内每一个键和键角，分别计算

### update.cpp
速度和位置更新，包括温度耦合的处理  
vv_update_1和2分别是velocity-verlet积分的第一步和第二步，原程序这里没有用蛙跳法，而是一个一阶精度算法，具体原因不明。力的计算应该夹在12两步中间，第二步的时候要不要更新加速度存疑。
### t_couple.cpp
温度耦合，先归约再求和，暂时不考虑归约的误差直接求和
不适配三维以外情况    
1、ekin_sum() 计算各个温度组的动能，得到半步温度Th和温度T，数组下标对应各个组  
    其中ekin_tot分了三个维度，是专门为后面压力耦合加上的
2、Berendsen_lambda() 由温度和参考温度ref_t计算温度耦合的参数lambda  
### virial.cpp
维里相关的参数和函数，主要用于压力耦合
为什么不直接与p_couple合并？区别在于这个要被include在各种力计算里，而p_couple相对独立一些
第一个版本里不考虑归约中的精度损失，直接逐个求和，那么二维数组是不必要的，直接用一个三元组即可
### vcm.cpp
消除质心运动，具体做的是对组内所有原子的动量求和，再除以总质量求得质心速度
每个group占据一个t_vcm结构体即可

### pretreat.cpp
前处理；输入力场文件、条件列表和结构文件，输出详细的可以直接计算的拓扑文件
csv读写系列函数。structure.csv对应gro文件。但gro文件也能带速度的，不确定我们有没有这个需求，暂时不实现
gro2csv: 很简单，就是gro文件转csv文件
t_atom_map：原子参数字典。读取%mole.itp和ffnonbonded.itp中的数据，构建参数字典  
            字典包含原子名称、原子序号、分子名称、原子质量、原子电荷、LJ参数
t_mole_map：分子参数字典。读取itp，给molecule结构分配成键、键角等参数
read_functype(): 读dump文件里的functype段  
    lj_para 2*n^2数组，按顺序存放c6和c12  
    ub_para 4*m 数组，按顺序存放theta,kTheta,r13,kUB四个参数  
### config.h
根据硬件架构不同使用的不同程序参数放在此文件中，自行调整
### vdwcf.cpp
范德华力和静电力  
float lj1 = lj_para[type];  
float lj2 = lj_para[type + 1];  
type均设定为偶数，前处理时从力场文件里查询两值  

### lincs.cpp
一般性约束算法
t_lincs_data:  
cons_num 约束的总数量  
A/B list 序号为约束序号，内容为键两端的原子序号  
len 序号为约束序号，内容为键长
Sdiag 序号为约束序号，内容为质量调和平均值 1 / sqrt （invmass atom1 + invmass atom2）  
dir_x/y/z 对应文献中的B矩阵，序号为约束序号，为键方向的单位向量  
nbr_cons 压缩矩阵，序号无意义，内容为某个约束的相邻约束  
coef 压缩矩阵，序号与上面的一一对应，atom_statics.m_inv[b] * lincs_data.Sdiag[i] * lincs_data.Sdiag[j];  
blcc 对应文献中的A矩阵，是约束间相互作用(constraint coupling),序号与上面的一一对应
nbr_cons_loc 上述压缩矩阵的索引，序号为约束序号(总容量为cons_num+1)


## 输入文件规则说明 
结构文件structure.csv  
例：  
分子序号 原子序号 分子名称 原子类型 x y z  
1,1,he,he,0.1,0.1,0.1  
2,1,he,he,0.2,0.2,0.2  
3,1,he,he,0.3,0.3,0.3  
与gro的主要区别在于全部用逗号隔开，不做对齐，容量较大  
方便用电子表格软件编辑  

