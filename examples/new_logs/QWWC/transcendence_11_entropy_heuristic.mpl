infolevel[Groebner]:=10;
Et_hat := [462532071134759-x_0, a*x_0-a*y_0-y_0*z_0+x_1, 27654837592025096323752656908-x_1, a*x_1+x_2+(-a-z_0)*y_1-y_0*z_1, -b*x_0-b*y_0+x_0*z_0+y_1, c*z_0+d*w_0-x_0*y_0+z_1, 2327261696006686982350191626670640107971252-x_2, a*x_2-2*z_1*y_1+x_3+(-a-z_0)*y_2-y_0*z_2, (-b+z_0)*x_1-b*y_1+y_2+x_0*z_1, c*z_1+d*w_1-x_0*y_1-x_1*y_0+z_2, -e*z_0+f*w_0-x_0*y_0+w_1, -15460657051372438519941789132700169821392928070012035129544-x_3, a*x_3-3*z_2*y_1-3*z_1*y_2+x_4+(-a-z_0)*y_3-y_0*z_3, 2*z_1*x_1+(-b+z_0)*x_2-b*y_2+y_3+x_0*z_2, c*z_2+d*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+z_3, -e*z_1+f*w_1-x_0*y_1-x_1*y_0+w_2, 2087310144461649513353687106909193349471468297179234169999889745921473824-x_4, a*x_4-4*z_3*y_1-6*z_2*y_2-4*z_1*y_3+x_5+(-a-z_0)*y_4-y_0*z_4, 3*z_2*x_1+3*z_1*x_2+(-b+z_0)*x_3-b*y_3+y_4+x_0*z_3, c*z_3+d*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+z_4, -e*z_2+f*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+w_3, 2898056534428371174690759898389123092528556169128723314534794667888624865342429389110704-x_5, a*x_5-5*z_4*y_1-10*z_3*y_2-10*z_2*y_3-5*z_1*y_4+x_6+(-a-z_0)*y_5-y_0*z_5, 4*z_3*x_1+6*z_2*x_2+4*z_1*x_3+(-b+z_0)*x_4-b*y_4+y_5+x_0*z_4, c*z_4+d*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+z_5, -e*z_3+f*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+w_4, 3311425124575611038181839911699882897574411313620363805054609360338040927269977340538148274079489124608-x_6, a*x_6-6*z_5*y_1-15*z_4*y_2-20*z_3*y_3-15*z_2*y_4-6*z_1*y_5+x_7+(-a-z_0)*y_6-y_0*z_6, 5*z_4*x_1+10*z_3*x_2+10*z_2*x_3+5*z_1*x_4+(-b+z_0)*x_5-b*y_5+y_6+x_0*z_5, c*z_5+d*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+z_6, -e*z_4+f*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+w_5, -64656120170310174326335330666443918894044955670610033436970907692423206666258657427341617525541510790555931549121504-x_7, a*x_7-7*z_6*y_1-21*z_5*y_2-35*z_4*y_3-35*z_3*y_4-21*z_2*y_5-7*z_1*y_6+x_8+(-a-z_0)*y_7-y_0*z_7, 6*z_5*x_1+15*z_4*x_2+20*z_3*x_3+15*z_2*x_4+6*z_1*x_5+(-b+z_0)*x_6-b*y_6+y_7+x_0*z_6, c*z_6+d*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+z_7, -e*z_5+f*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+w_6, -9908573056011909293045330416744625761402810808186732533394669548537205798315191294901472565828263434059357895268617610447376269821088-x_8, a*x_8-8*z_7*y_1-28*z_6*y_2-56*z_5*y_3-70*z_4*y_4-56*z_3*y_5-28*z_2*y_6-8*z_1*y_7+x_9+(-a-z_0)*y_8-y_0*z_8, 7*z_6*x_1+21*z_5*x_2+35*z_4*x_3+35*z_3*x_4+21*z_2*x_5+7*z_1*x_6+(-b+z_0)*x_7-b*y_7+y_8+x_0*z_7, c*z_7+d*w_7-x_0*y_7-7*x_1*y_6-21*x_2*y_5-35*x_3*y_4-35*x_4*y_3-21*x_5*y_2-7*x_6*y_1-x_7*y_0+z_8, -e*z_6+f*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+w_7, -392286824968509616140136077716056452945139281325461800474719899044269163701254556797938905269459610506485936237428706991873142274092837429448057600-x_9, 6349572714173871959973839081256977289551260116401318207620897387440769383089668654070273677425207455798-w_6, z_aux-1];
vars:=[x_9, z_8, y_8, x_8, z_7, y_7, x_7, w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, z_4, y_4, x_4, w_4, z_3, y_3, x_3, w_3, z_2, y_2, x_2, w_2, z_1, y_1, x_1, w_1, z_0, y_0, x_0, w_0, z_aux, w_aux, a, b, c, d, e, f];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [w_6] 5
# 228, [1 = 8, -x_4*y_0 = 2, -x_0*y_2 = 2, -x_0*y_6 = 2, -2*x_1*y_1 = 2, -x_5*y_0 = 2, -x_0*y_3 = 2, -6*x_1*y_5 = 2, -x_0*y_1 = 2, -x_0*y_0 = 2, -3*x_1*y_2 = 2, -4*x_1*y_3 = 2, -15*x_2*y_4 = 2, -3*x_2*y_1 = 2, -4*x_3*y_1 = 2, -20*x_3*y_3 = 2, -x_1*y_0 = 2, -6*x_2*y_2 = 2, -x_3*y_0 = 2, -5*x_1*y_4 = 2, -15*x_4*y_2 = 2, -x_6*y_0 = 2, -x_0*y_4 = 2, -6*x_5*y_1 = 2, -x_0*y_5 = 2, -5*x_4*y_1 = 2, -x_2*y_0 = 2, -10*x_2*y_3 = 2, -10*x_3*y_2 = 2, -15*z_4*y_2 = 1, x_7 = 1, -b*y_5 = 1, -35*x_4*y_3 = 1, z_2*x_5 = 1, z_2*x_3 = 1, -x_7*b = 1, z_aux = 1, -y_0*z_0 = 1, y_1 = 1, w_2 = 1, -3*z_2*y_1 = 1, z_1 = 1, -x_2 = 1, -y_5*a = 1, -35*x_3*y_4 = 1, -138067306532032552105501511594993580674626965293010153513125713762985305968281696719375309636542612405007387700132638946 = 1, -6*z_2*y_2 = 1, z_5*x_2 = 1, -a*y_0 = 1, d*w_0 = 1, z_7 = 1, -e*z_4 = 1, -x_5*b = 1, -b*y_7 = 1, -x_7 = 1, -b*y_1 = 1, a*x_2 = 1, -x_5 = 1, -y_0*z_4 = 1, x_2 = 1, z_4*x_2 = 1, z_3*x_3 = 1, x_9 = 1, -6*z_5*y_1 = 1, -y_0*z_3 = 1, y_3 = 1, z_4 = 1, w_1 = 1, -y_0*z_6 = 1, -y_7*z_0 = 1, -e*z_3 = 1, -x_1*b = 1, x_0*z_4 = 1, y_4 = 1, z_3 = 1, -21*x_2*y_5 = 1, -4*z_3*y_1 = 1, f*w_5 = 1, -y_4*a = 1, z_1*x_1 = 1, z_1*x_5 = 1, -x_2*b = 1, z_1*x_6 = 1, f*w_2 = 1, -5*z_4*y_1 = 1, x_5 = 1, -7*x_1*y_6 = 1, -21*z_2*y_5 = 1, -8*z_1*y_7 = 1, a*x_1 = 1, z_8 = 1, -y_6*a = 1, f*w_4 = 1, -7*z_1*y_6 = 1, c*z_3 = 1, x_4*z_0 = 1, x_1 = 1, x_2*z_0 = 1, z_2*x_2 = 1, w_5 = 1, a*x_8 = 1, c*z_2 = 1, -x_0*y_7 = 1, -35*z_3*y_4 = 1, c*z_1 = 1, -b*y_0 = 1, x_7*z_0 = 1, a*x_7 = 1, a*x_5 = 1, -b*y_2 = 1, a*x_6 = 1, c*z_4 = 1, x_0*z_7 = 1, x_3*z_0 = 1, -y_3*z_0 = 1, f*w_3 = 1, f*w_1 = 1, -x_0 = 1, -y_8*z_0 = 1, z_6*x_1 = 1, w_4 = 1, -3*z_1*y_2 = 1, -y_0*z_5 = 1, -b*x_0 = 1, d*w_7 = 1, z_2 = 1, -20*z_3*y_3 = 1, -21041093569179189982427222998 = 1, f = 1, -y_1*z_0 = 1, -5*z_1*y_4 = 1, -70*z_4*y_4 = 1, -21*z_5*y_2 = 1, -e*z_2 = 1, d*w_1 = 1, -7*z_6*y_1 = 1, d*w_2 = 1, x_5*z_0 = 1, -b*y_3 = 1, -x_9 = 1, -x_6 = 1, x_6*z_0 = 1, -33105577711812079681034496076987347178833991680801882791801999330087088606 = 1, -21*x_5*y_2 = 1, z_2*x_4 = 1, -y_2*a = 1, z_4*x_3 = 1, z_1*x_3 = 1, -7*x_6*y_1 = 1, -56*z_3*y_5 = 1, -35*z_4*y_3 = 1, z_6 = 1, y_2 = 1, x_4 = 1, y_7 = 1, -10*z_2*y_3 = 1, -x_3*b = 1, -15*z_2*y_4 = 1, c*z_7 = 1, -x_4 = 1, -e*z_0 = 1, w_3 = 1, -x_8 = 1, -6*z_1*y_5 = 1, z_2*x_1 = 1, -y_0*z_1 = 1, y_8 = 1, -28*z_2*y_6 = 1, x_1*z_0 = 1, z_3*x_4 = 1, -10*z_3*y_2 = 1, z_4*x_1 = 1, c*z_6 = 1, c*z_5 = 1, -y_0*z_2 = 1, -e*z_1 = 1, -x_7*y_0 = 1, -e*z_5 = 1, a*x_0 = 1, w_7 = 1, -x_4*b = 1, c*z_0 = 1, d = 1, -8*z_7*y_1 = 1, -y_7*a = 1, z_5 = 1, z_1*x_2 = 1, -y_1*a = 1, -28*z_6*y_2 = 1, z_5*x_1 = 1, -y_6*z_0 = 1, -y_4*z_0 = 1, x_0*z_1 = 1, -y_3*a = 1, -e*z_6 = 1, a*x_3 = 1, -2*z_1*y_1 = 1, z_3*x_2 = 1, -b*y_4 = 1, x_0*z_6 = 1, -x_1 = 1, -b*y_6 = 1, a*x_4 = 1, x_0*z_2 = 1, d*w_5 = 1, x_0*z_5 = 1, -x_3 = 1, -y_0*z_8 = 1, z_3*x_1 = 1, -y_0*z_7 = 1, y_6 = 1, x_0*z_3 = 1, y_5 = 1, -y_2*z_0 = 1, -4*z_1*y_3 = 1, -y_8*a = 1, x_8 = 1, d*w_3 = 1, z_1*x_4 = 1, x_3 = 1, -56*z_5*y_3 = 1, x_0*z_0 = 1, d*w_4 = 1, -y_5*z_0 = 1, -x_6*b = 1, f*w_0 = 1, -1 = 1, x_6 = 1]
# 263, -5.361310784
# 29
# [w_6 = 3]
# [w_6 = [2, 1, 2]]