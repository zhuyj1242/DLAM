import numpy as np
fn = 16
# 层数
ln = 16
# 半跨段数+1
h = 160
L = 2550
Ec = 34000
Es = 200000
n, m, k = 1, 0, 0
Kt = 2.2e9
fm = 0.00001
fc = 50
As1, As2 = 633.6, 633.6
et = -np.ones((fn, ln))/10000
expand = 3.7
aerfa = np.ones(ln)
disp1 = 80
# 5,10,15,20,30,40,50,60,70,
K1 = 140000
loop_disp = 0
wi = []
wj = []
Tn = []
Tnerror = []
expanerr = []
Merror = []
werror = []
Tnn = []
Load = []
for disp in disp1:
    loop_disp = loop_disp + 1
    loop_k = 0
    for K in K1:
        loop_k = loop_k + 1
        upboundf = 0.0006
        lowboundf = 0
        loopfm = 0
        while 1
            loopfm = loopfm + 1
            upboundaerfa= np.ones((1, ln))
            lowboundaerfa=-np.ones((1, ln))
            loopaerfa = 0
            if loopfm == 80
                break
            while 1
                loopaerfa = loopaerfa + 1
                for j in range(ln):
                    fi[j] = aerfa[j] * fm
                    upboundet = 0.002
                    lowboundet = -0.002
                    loopet = 0
                    while 1
                        if loopet == 50:
                            break
                        loopet = loopet + 1
                        Tcn[j] = 0
                        for i in range(1,fn):
                            et[i, j] = et[1, j] + fi[j] * h / fn * [i - 1] # 正弯矩曲率是正值
                        for i in range(fn):
                            if et(i, j) > 0:
                                sigmac[i, j] = 0
                            elif - 0.0004 < et[i, j] < 0:
                                sigmac[i, j] = et[i, j] * Ec
                            elif - 0.002 < et[i, j] < -0.0004:
                                sigmac[i, j] = -(2 * et[i, j] / (-0.002) - (et[i, j] / (-0.002)) ^ 2) * fc
                            else:
                                sigmac[i, j] = et[i, j] * 0.8 * Ec
                            Tcn[j] = Tcn[j] + sigmac[i, j] * h / fn * 1000
                        if et[fn, j] > 0.002:
                            sigmas[j] = 0.002 * Es
                        elif et[fn, j] < -0.002:
                            sigmas[j] = -0.002 * Es
                        else:
                            sigmas[j] = et[fn, j] * Es
                        if et[1, j] > 0.002:
                            sigmass[j] = 0.002 * Es
                        elif et[1, j] < -0.002:
                            sigmass[j] = -0.002 * Es
                        else:
                            sigmass[j] = et[1, j] * Es
                        Ts[j] = sigmas[j] * As1
                        Tss[j] = sigmass[j] * As2
                        Tn[j] = Tcn[j] + Ts[j] + Tss[j] # 截面总轴力
                        Tn1 = -K * expand # 这里加了个负号，轴力为负
                        if abs(Tn[j] - Tn1) < abs(Tn1) / 50:
                            break
                        elif Tn[j] - Tn1 > 0:
                            upboundet = et[1, j]
                            et[1, j] = (upboundet + lowboundet) / 2
                        else:
                            lowboundet = et[1, j]
                            et[1, j] = (upboundet + lowboundet) / 2
# 完成截面应变假设
                    Mc = []
                    for i in range (fn):
                        Mc[j] = Mc[j] - sigmac[i, j] * h / fn * (h / 2 - h / fn * (i - 0.5)) * 1000 # 混凝土部分的弯矩
                    M[j] = Mc[j] + Ts[j] * h * (0.5 - 1 / fn / 2) - Tss[j] * h * (0.5 - 1 / fn / 2) # 假设钢筋位于最边缘层的中心？
# 计算梁端转角
                angle = 0
                for jj in range (1, ln):
                    angle = angle + (fi[jj - 1] + fi[jj]) / 2 * L / 2 / (ln - 1)
# 计算挠度，正值表示下挠
                wi[ln] = 0
                for jj in range (1, ln):
                    wi[ln] = wi[ln] + (fi[jj - 1] + fi[jj]) / 2 * L / 2 / (ln - 1) * ((jj - 1.5) * L / 2 / (ln - 1)) # 这个是对了，论文里和这个是不是一致再检查一下，可以把坐标也换过来
                for ii in range (ln-1):
                    wj[ii] = 0
                    for jj in range (ii, ln):
                        wj[ii] = wj[ii] + (fi[jj - 1] + fi[jj]) / 2 * L / 2 / (ln - 1) * ((jj - ii - 0.5) * L / 2 / (ln - 1)) # 好绕啊，这个也是对的
                    wi[ii] = wi[ln] - wj[ii]
# 基于梁端转角和挠度，计算弯矩分布
                for i in range (ln):
                    Mi[i] = -Kt * angle + (M[ln] + K * expand * (h / 2 - wi[ln]) + Kt * angle) * ((i - 1) / (ln - 1)) - K * expand * (h / 2 - wi[i]) # 和论文里符号不一致，这里将来要按论文里换成eL，目前是按EL全长均为0.5h计算的
# 基于弯矩分布重调曲率
                q = 0
                for i in range (ln-1):
                    if abs(M[i] - Mi[i]) < abs(Mi[i]) / 50:
                        q = q + 1 # q来判断弯矩满足要求的个数
                        continue
                    elif M[i] - Mi[i] > 0:
                        upboundaerfa[i] = aerfa[i]
                        aerfa[i] = (upboundaerfa[i] + lowboundaerfa[i]) / 2
                    else:
                        lowboundaerfa[i] = aerfa[i]
                        aerfa[i] = (upboundaerfa[i] + lowboundaerfa[i]) / 2
                if q == ln - 1:
                    break
                if loopaerfa == 80:
                    break
# 进行跨中曲率修正
            if abs(wi[ln] - disp) < (disp / 20):
                break
            elif wi[ln] - disp > 0:
                upboundf = fm
                fm = (upboundf + lowboundf) / 2
            else:
                lowboundf = fm
                fm = (upboundf + lowboundf) / 2
# 修正压杆力
        expand1 = 0
        for i in range (1,ln):
            expand1 = expand1 + (et[fn, i - 1] + et[fn, i]) / 2 * L / 2 / (ln - 1)

        Tnerror[loop_disp,:]=(Tn - Tn1) / Tn1 * 100
        expanerr[loop_disp] = (expand1 - expand) / expand
        Merror[loop_disp,:]=(Mi - M) / M[ln] * 1000
        werror[loop_disp] = (wi[ln] - disp) / disp
        Tnn[loop_disp, loop_k] = Tn1
        Load[loop_disp, loop_k] = 4 / L * (M[ln] + K * expand * (h / 2 - wi[ln]) + Kt * angle)

#figure;
#plot(disp1, Load(:, 1) / 1000, 's-', 'MarkerSize', 3, 'Linewidth', 0.7)
# hold
#on;
# plot(disp1, Load(:, 2), 'o-', 'MarkerSize', 3, 'Linewidth', 0.7)
# hold
#on;
# plot(disp1, Load(:, 3), 'o-', 'MarkerSize', 3, 'Linewidth', 0.7)
#xlabel('disp(mm)');
#ylabel('load(kN)');
# h1 = legend('Test', 'CEP90', 'location', 'southeast'); % 'ACI', 'GL2000', 'B4s',
# set(h1, 'edgecolor', 'none');

#set(gcf, 'Units', 'centimeters', 'Position', [7 8 6.8 5.6]);
#set(gca, 'ticklength', [0.011 0.05], 'Linewidth', 0.7, 'fontsize', 9, 'fontname', 'Times New Roman'); %, 'ylim', [0, 3]
#set(gca, 'ticklength', [0.011 0.05], 'Linewidth', 0.7, 'fontsize', 9, 'fontname', 'Times New Roman');

