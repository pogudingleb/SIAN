infolevel[Groebner]:=10;
Et_hat := [60015158889925-x_0, a*x_0-a*y_0-y_0*z_0+x_1, 482794927070516-e_0, e_1, 435168281105208313529470780820-x_1, a*x_1+x_2+(-a-z_0)*y_1-y_0*z_1, -b*x_0-b*y_0+x_0*z_0+y_1, c*z_0+d*w_0-x_0*y_0+z_1, -349056515701663558208448694275409940338695145-x_2, a*x_2-2*z_1*y_1+x_3+(-a-z_0)*y_2-y_0*z_2, (-b+z_0)*x_1-b*y_1+y_2+x_0*z_1, c*z_1+d*w_1-x_0*y_1-x_1*y_0+z_2, -e_0*z_0+f*w_0-x_0*y_0+w_1, 342580187474050372148232963014722671116145534925971935916965-x_3, a*x_3-3*z_2*y_1-3*z_1*y_2+x_4+(-a-z_0)*y_3-y_0*z_3, 2*z_1*x_1+(-b+z_0)*x_2-b*y_2+y_3+x_0*z_2, c*z_2+d*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+z_3, -e_0*z_1-e_1*z_0+f*w_1-x_0*y_1-x_1*y_0+w_2, 894946964074237983602597614808964899598723784267796004645251574996792444680-x_4, a*x_4-4*z_3*y_1-6*z_2*y_2-4*z_1*y_3+x_5+(-a-z_0)*y_4-y_0*z_4, 3*z_2*x_1+3*z_1*x_2+(-b+z_0)*x_3-b*y_3+y_4+x_0*z_3, c*z_3+d*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+z_4, -e_0*z_2-2*e_1*z_1-e_2*z_0+f*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+w_3, e_2, -6179585922441252272915657701908999267891853241947802751952052646525676203218621640743154985-x_5, a*x_5-5*z_4*y_1-10*z_3*y_2-10*z_2*y_3-5*z_1*y_4+x_6+(-a-z_0)*y_5-y_0*z_5, 4*z_3*x_1+6*z_2*x_2+4*z_1*x_3+(-b+z_0)*x_4-b*y_4+y_5+x_0*z_4, c*z_4+d*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+z_5, -e_0*z_3-3*e_1*z_2-3*e_2*z_1-e_3*z_0+f*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+w_4, e_3, 25184325322525392347137611579479668617946817224535813775252355481999903718867820974297123117545122670020905-x_6, a*x_6-6*z_5*y_1-15*z_4*y_2-20*z_3*y_3-15*z_2*y_4-6*z_1*y_5+x_7+(-a-z_0)*y_6-y_0*z_6, 5*z_4*x_1+10*z_3*x_2+10*z_2*x_3+5*z_1*x_4+(-b+z_0)*x_5-b*y_5+y_6+x_0*z_5, c*z_5+d*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+z_6, -e_0*z_4-4*e_1*z_3-6*e_2*z_2-4*e_3*z_1-e_4*z_0+f*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+w_5, e_4, -50708780792545781707902930423754438518570974216643656411434000103872558592248841297907190822258118414078177554593613957700-x_7, a*x_7-7*z_6*y_1-21*z_5*y_2-35*z_4*y_3-35*z_3*y_4-21*z_2*y_5-7*z_1*y_6+x_8+(-a-z_0)*y_7-y_0*z_7, 6*z_5*x_1+15*z_4*x_2+20*z_3*x_3+15*z_2*x_4+6*z_1*x_5+(-b+z_0)*x_6-b*y_6+y_7+x_0*z_6, c*z_6+d*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+z_7, -e_0*z_5-5*e_1*z_4-10*e_2*z_3-10*e_3*z_2-5*e_4*z_1-e_5*z_0+f*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+w_6, e_5, -21134025251868778549702022939993186706611735389572075068855207936537685213313103515845423301101867204814634416762355051237424514853177525-x_8, a*x_8-8*z_7*y_1-28*z_6*y_2-56*z_5*y_3-70*z_4*y_4-56*z_3*y_5-28*z_2*y_6-8*z_1*y_7+x_9+(-a-z_0)*y_8-y_0*z_8, 7*z_6*x_1+21*z_5*x_2+35*z_4*x_3+35*z_3*x_4+21*z_2*x_5+7*z_1*x_6+(-b+z_0)*x_7-b*y_7+y_8+x_0*z_7, c*z_7+d*w_7-x_0*y_7-7*x_1*y_6-21*x_2*y_5-35*x_3*y_4-35*x_4*y_3-21*x_5*y_2-7*x_6*y_1-x_7*y_0+z_8, -e_0*z_6-6*e_1*z_5-15*e_2*z_4-20*e_3*z_3-15*e_4*z_2-6*e_5*z_1-e_6*z_0+f*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+w_7, e_6, 1129886440817026226645242644370512229284744464268923555587867155509021139283795180699128724322442061350018351952005093981043740158654845005157050210530725-x_9, -e_1, -e_2, -e_3, -e_4, -e_5, -e_6, z_aux-1];
vars:=[x_9, z_8, y_8, x_8, z_7, y_7, x_7, w_7, z_6, y_6, x_6, w_6, e_6, z_5, y_5, x_5, w_5, e_5, z_4, y_4, x_4, w_4, e_4, z_3, y_3, x_3, w_3, e_3, z_2, y_2, x_2, w_2, e_2, z_1, y_1, x_1, w_1, e_1, z_0, y_0, x_0, w_0, e_0, z_aux, w_aux, a, b, c, d, f];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [e] 14
# 226, [1 = 10, -6*x_1*y_5 = 2, -4*x_1*y_3 = 2, -5*x_1*y_4 = 2, -x_0*y_5 = 2, -x_0*y_6 = 2, -x_0*y_4 = 2, -10*x_2*y_3 = 2, -15*x_2*y_4 = 2, -10*x_3*y_2 = 2, -x_5*y_0 = 2, -5*x_4*y_1 = 2, -3*x_1*y_2 = 2, -x_2*y_0 = 2, -x_0*y_1 = 2, -x_0*y_3 = 2, -x_6*y_0 = 2, -2*x_1*y_1 = 2, -x_4*y_0 = 2, -x_0*y_2 = 2, -20*x_3*y_3 = 2, -4*x_3*y_1 = 2, -x_3*y_0 = 2, -15*x_4*y_2 = 2, -x_0*y_0 = 2, -3*x_2*y_1 = 2, -x_1*y_0 = 2, -6*x_2*y_2 = 2, -6*x_5*y_1 = 2, -y_1*z_0 = 1, d*w_2 = 1, a*x_2 = 1, -x_7*b = 1, d*w_3 = 1, -2*z_1*y_1 = 1, -x_4 = 1, -z_0 = 1, c*z_7 = 1, -7*z_1*y_6 = 1, z_1*x_3 = 1, -10*z_2*y_3 = 1, -x_3*b = 1, w_3 = 1, -35*x_4*y_3 = 1, z_4*x_1 = 1, x_7 = 1, -b*y_3 = 1, c*z_6 = 1, -21*z_2*y_5 = 1, -21*x_5*y_2 = 1, -b*y_7 = 1, z_2 = 1, z_2*x_1 = 1, -y_0*z_1 = 1, x_0*z_1 = 1, c*z_5 = 1, -6*z_1*y_5 = 1, a*x_6 = 1, f*w_4 = 1, -15*z_2*y_4 = 1, -x_9 = 1, -z_3 = 1, -z_5 = 1, -5*z_1*y_4 = 1, -y_0*z_2 = 1, x_5 = 1, -x_3 = 1, x_6*z_0 = 1, -b*x_0 = 1, x_4 = 1, d*w_5 = 1, z_1*x_1 = 1, x_1*z_0 = 1, -y_2*a = 1, a*x_8 = 1, -x_6 = 1, y_2 = 1, -7*x_6*y_1 = 1, z_2*x_5 = 1, -z_4 = 1, x_1 = 1, -y_5*z_0 = 1, -35*z_3*y_4 = 1, c*z_1 = 1, -x_1*b = 1, z_1*x_4 = 1, -z_1 = 1, z_2*x_3 = 1, -x_5*b = 1, a*x_4 = 1, z_4*x_3 = 1, -y_4*z_0 = 1, z_3*x_3 = 1, -x_5 = 1, -21*x_2*y_5 = 1, -x_7 = 1, x_0*z_5 = 1, y_3 = 1, d*w_6 = 1, a*x_7 = 1, y_4 = 1, x_3*z_0 = 1, -y_7*a = 1, y_5 = 1, z_3*x_1 = 1, -x_4*b = 1, d*w_1 = 1, -35*x_3*y_4 = 1, -56*z_5*y_3 = 1, x_3 = 1, d*w_7 = 1, -y_7*z_0 = 1, x_0*z_4 = 1, -70*z_4*y_4 = 1, c*z_3 = 1, -y_2*z_0 = 1, -b*y_4 = 1, f*w_5 = 1, d*w_4 = 1, y_7 = 1, -y_3*a = 1, a*x_1 = 1, d*w_0 = 1, -y_6*z_0 = 1, -y_1*a = 1, -6*z_5*y_1 = 1, w_2 = 1, -y_0*z_3 = 1, a*x_3 = 1, z_2*x_2 = 1, z_5*x_2 = 1, z_2*x_4 = 1, -z_6 = 1, f*w_0 = 1, c*z_2 = 1, y_1 = 1, w_5 = 1, -4*z_3*y_1 = 1, a*x_0 = 1, -x_1 = 1, x_0*z_6 = 1, -y_0*z_7 = 1, -3*z_2*y_1 = 1, -y_0*z_8 = 1, x_0*z_2 = 1, -x_2 = 1, c*z_0 = 1, f*w_6 = 1, -8*z_7*y_1 = 1, -7*x_1*y_6 = 1, -y_0*z_0 = 1, -4*z_1*y_3 = 1, -56*z_3*y_5 = 1, x_0*z_3 = 1, w_7 = 1, x_9 = 1, -x_6*b = 1, z_6 = 1, -7*z_6*y_1 = 1, -b*y_5 = 1, z_8 = 1, z_3*x_4 = 1, z_4*x_2 = 1, -15*z_4*y_2 = 1, -b*y_1 = 1, -28*z_6*y_2 = 1, -y_8*a = 1, -20*z_3*y_3 = 1, z_1*x_6 = 1, -3*z_1*y_2 = 1, x_2*z_0 = 1, z_7 = 1, -x_0*y_7 = 1, -b*y_6 = 1, f*w_2 = 1, c*z_4 = 1, -y_5*a = 1, f*w_3 = 1, -6*z_2*y_2 = 1, -5*z_4*y_1 = 1, -z_2 = 1, -y_3*z_0 = 1, x_0*z_0 = 1, x_0*z_7 = 1, x_8 = 1, z_3 = 1, w_1 = 1, -x_7*y_0 = 1, -y_4*a = 1, -y_8*z_0 = 1, -x_8 = 1, z_1*x_5 = 1, -b*y_0 = 1, x_7*z_0 = 1, z_1*x_2 = 1, z_5 = 1, z_5*x_1 = 1, y_8 = 1, -21*z_5*y_2 = 1, x_2 = 1, -a*y_0 = 1, -y_6*a = 1, y_6 = 1, -y_0*z_5 = 1, x_6 = 1, z_aux = 1, f*w_1 = 1, -b*y_2 = 1, z_1 = 1, -8*z_1*y_7 = 1, -y_0*z_4 = 1, w_4 = 1, z_4 = 1, a*x_5 = 1, -10*z_3*y_2 = 1, -x_2*b = 1, -y_0*z_6 = 1, z_6*x_1 = 1, x_4*z_0 = 1, -x_0 = 1, -35*z_4*y_3 = 1, x_5*z_0 = 1, z_3*x_2 = 1, -1 = 1, -28*z_2*y_6 = 1, w_6 = 1]
# 263, -5.337012995
# 29
# [e = 7]
# [e = [2, 2, 2, 2, 2, 2, 2]]