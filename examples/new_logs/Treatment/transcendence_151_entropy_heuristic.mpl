infolevel[Groebner]:=10;
Et_hat := [35313285411072-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 173491112598457-b_0, b_1, 109231573633617-d_0, d_1, -Tr_1+22429201749942712379271892704, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b_0*d_0+In_0*N*a+In_0*N*g-In_0*S_0*b_0+In_1*N, -Tr_2+180286390458581285282347417207295652725588853730283328292329620375712/53719254054335, -In_2*g+Tr_2*nu+Tr_3, ((-Tr_0*d_1-Tr_1*d_0-In_1)*b_0-b_1*(Tr_0*d_0+In_0))*S_0-S_1*(Tr_0*d_0+In_0)*b_0+((a+g)*In_1+In_2)*N, S_0*Tr_0*b_0*d_0+In_0*S_0*b_0+N*S_1, -Tr_3-8043317859371202393416290788006257758627023368199437268695928336391482358110393498609899937840581463607078752/577151651230837467026458445, -In_3*g+Tr_3*nu+Tr_4, ((-Tr_0*d_2-2*Tr_1*d_1-Tr_2*d_0-In_2)*b_0+(-2*b_1*d_1-b_2*d_0)*Tr_0-2*Tr_1*b_1*d_0-In_0*b_2-2*In_1*b_1)*S_0+((-2*S_1*d_1-S_2*d_0)*Tr_0-2*d_0*Tr_1*S_1-S_2*In_0-2*S_1*In_1)*b_0-2*S_1*Tr_0*b_1*d_0-2*In_0*S_1*b_1+N*((a+g)*In_2+In_3), ((Tr_0*d_1+Tr_1*d_0+In_1)*b_0+b_1*(Tr_0*d_0+In_0))*S_0+S_1*(Tr_0*d_0+In_0)*b_0+N*S_2, b_2, d_2, -Tr_4+8971137813310425175997330047331564746121436640989991620493489351543446446677946126320683335974828637249702626514041165217741530896583716971087320698144/155020780901741527467692498491413248045375, -In_4*g+Tr_4*nu+Tr_5, ((-Tr_0*d_3-3*Tr_1*d_2-3*Tr_2*d_1-Tr_3*d_0-In_3)*b_0+(-3*b_1*d_2-3*b_2*d_1-b_3*d_0)*Tr_0+(-3*Tr_1*b_2-3*Tr_2*b_1)*d_0+(-6*Tr_1*d_1-3*In_2)*b_1-3*In_1*b_2-b_3*In_0)*S_0+((-3*S_1*d_2-3*S_2*d_1-S_3*d_0)*Tr_0+(-3*S_1*Tr_2-3*S_2*Tr_1)*d_0+(-6*Tr_1*d_1-3*In_2)*S_1-S_3*In_0-3*In_1*S_2)*b_0+((-3*S_1*b_2-3*S_2*b_1)*d_0-6*d_1*S_1*b_1)*Tr_0-6*S_1*Tr_1*b_1*d_0+(-3*In_0*b_2-6*In_1*b_1)*S_1-3*In_0*S_2*b_1+((a+g)*In_3+In_4)*N, ((Tr_0*d_2+2*Tr_1*d_1+Tr_2*d_0+In_2)*b_0+(2*b_1*d_1+b_2*d_0)*Tr_0+2*Tr_1*b_1*d_0+In_0*b_2+2*In_1*b_1)*S_0+((2*S_1*d_1+S_2*d_0)*Tr_0+2*d_0*Tr_1*S_1+S_2*In_0+2*S_1*In_1)*b_0+2*S_1*Tr_0*b_1*d_0+2*In_0*S_1*b_1+N*S_3, b_3, d_3, -Tr_5-400239379184194246327767248369432440048994248264646030465679028036786434679640942982328216351995186132574506860640096448320104672726141587938249203785435058282537283026355923975521705665801312/1665520142592411257336215203062383399249088396559090125, -In_5*g+Tr_5*nu+Tr_6, ((-Tr_0*d_4-4*Tr_1*d_3-6*Tr_2*d_2-4*Tr_3*d_1-Tr_4*d_0-In_4)*b_0+(-4*b_1*d_3-6*b_2*d_2-4*b_3*d_1-b_4*d_0)*Tr_0+(-4*Tr_1*b_3-6*Tr_2*b_2-4*Tr_3*b_1)*d_0+(-12*Tr_1*d_2-12*Tr_2*d_1-4*In_3)*b_1-12*b_2*d_1*Tr_1-6*b_2*In_2-4*b_3*In_1-b_4*In_0)*S_0+((-4*S_1*d_3-6*S_2*d_2-4*S_3*d_1-S_4*d_0)*Tr_0+(-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1)*d_0+(-12*Tr_1*d_2-12*Tr_2*d_1-4*In_3)*S_1-12*S_2*d_1*Tr_1-In_0*S_4-4*In_1*S_3-6*In_2*S_2)*b_0+((-4*S_1*b_3-6*S_2*b_2-4*S_3*b_1)*d_0+(-12*b_1*d_2-12*b_2*d_1)*S_1-12*d_1*S_2*b_1)*Tr_0+((-12*Tr_1*b_2-12*Tr_2*b_1)*S_1-12*Tr_1*S_2*b_1)*d_0+((-24*Tr_1*d_1-12*In_2)*b_1-12*In_1*b_2-4*b_3*In_0)*S_1+(-4*In_0*S_3-12*In_1*S_2)*b_1-6*In_0*S_2*b_2+((a+g)*In_4+In_5)*N, ((Tr_0*d_3+3*Tr_1*d_2+3*Tr_2*d_1+Tr_3*d_0+In_3)*b_0+(3*b_1*d_2+3*b_2*d_1+b_3*d_0)*Tr_0+(3*Tr_1*b_2+3*Tr_2*b_1)*d_0+(6*Tr_1*d_1+3*In_2)*b_1+b_3*In_0+3*In_1*b_2)*S_0+((3*S_1*d_2+3*S_2*d_1+S_3*d_0)*Tr_0+(3*S_1*Tr_2+3*S_2*Tr_1)*d_0+(6*Tr_1*d_1+3*In_2)*S_1+S_3*In_0+3*In_1*S_2)*b_0+((3*S_1*b_2+3*S_2*b_1)*d_0+6*d_1*S_1*b_1)*Tr_0+6*S_1*Tr_1*b_1*d_0+(3*In_0*b_2+6*In_1*b_1)*S_1+3*In_0*S_2*b_1+N*S_4, b_4, d_4, -Tr_6+446408147949901437226903545584052789333599711958020848319837133039514714028294660758826504962521448326552904483517973570017486641120615267847308260833739589837100684079565217172839348881308503895213440242845794299674664554311011667232/447352498362669978765310867796665023640198954207048532304024559709375, -In_6*g+Tr_6*nu+Tr_7, ((-Tr_0*d_5-5*Tr_1*d_4-10*Tr_2*d_3-10*Tr_3*d_2-5*Tr_4*d_1-Tr_5*d_0-In_5)*b_0+(-5*b_1*d_4-10*b_2*d_3-10*b_3*d_2-5*b_4*d_1-b_5*d_0)*Tr_0+(-5*Tr_1*b_4-10*Tr_2*b_3-10*Tr_3*b_2-5*Tr_4*b_1)*d_0+(-20*Tr_1*d_3-30*Tr_2*d_2-20*Tr_3*d_1-5*In_4)*b_1+(-30*b_2*d_2-20*b_3*d_1)*Tr_1-30*Tr_2*b_2*d_1-10*b_2*In_3-10*In_2*b_3-5*In_1*b_4-In_0*b_5)*S_0+((-5*S_1*d_4-10*S_2*d_3-10*S_3*d_2-5*S_4*d_1-S_5*d_0)*Tr_0+(-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1)*d_0+(-20*Tr_1*d_3-30*Tr_2*d_2-20*Tr_3*d_1-5*In_4)*S_1+(-30*S_2*d_2-20*S_3*d_1)*Tr_1-30*d_1*Tr_2*S_2-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2)*b_0+((-5*S_1*b_4-10*S_2*b_3-10*S_3*b_2-5*S_4*b_1)*d_0+(-20*b_1*d_3-30*b_2*d_2-20*b_3*d_1)*S_1+(-30*S_2*d_2-20*S_3*d_1)*b_1-30*b_2*d_1*S_2)*Tr_0+((-20*Tr_1*b_3-30*Tr_2*b_2-20*Tr_3*b_1)*S_1+(-30*S_2*Tr_2-20*S_3*Tr_1)*b_1-30*b_2*Tr_1*S_2)*d_0+((-60*Tr_1*d_2-60*Tr_2*d_1-20*In_3)*b_1-60*b_2*d_1*Tr_1-30*b_2*In_2-20*b_3*In_1-5*b_4*In_0)*S_1+(-60*S_2*Tr_1*d_1-5*In_0*S_4-20*In_1*S_3-30*In_2*S_2)*b_1+(-10*In_0*b_3-30*In_1*b_2)*S_2-10*In_0*S_3*b_2+((a+g)*In_5+In_6)*N, ((Tr_0*d_4+4*Tr_1*d_3+6*Tr_2*d_2+4*Tr_3*d_1+Tr_4*d_0+In_4)*b_0+(4*b_1*d_3+6*b_2*d_2+4*b_3*d_1+b_4*d_0)*Tr_0+(4*Tr_1*b_3+6*Tr_2*b_2+4*Tr_3*b_1)*d_0+(12*Tr_1*d_2+12*Tr_2*d_1+4*In_3)*b_1+12*b_2*d_1*Tr_1+6*b_2*In_2+4*b_3*In_1+b_4*In_0)*S_0+((4*S_1*d_3+6*S_2*d_2+4*S_3*d_1+S_4*d_0)*Tr_0+(4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1)*d_0+(12*Tr_1*d_2+12*Tr_2*d_1+4*In_3)*S_1+12*S_2*d_1*Tr_1+In_0*S_4+4*In_1*S_3+6*In_2*S_2)*b_0+((4*S_1*b_3+6*S_2*b_2+4*S_3*b_1)*d_0+(12*b_1*d_2+12*b_2*d_1)*S_1+12*d_1*S_2*b_1)*Tr_0+((12*Tr_1*b_2+12*Tr_2*b_1)*S_1+12*Tr_1*S_2*b_1)*d_0+((24*Tr_1*d_1+12*In_2)*b_1+12*In_1*b_2+4*b_3*In_0)*S_1+(4*In_0*S_3+12*In_1*S_2)*b_1+6*In_0*S_2*b_2+N*S_5, b_5, d_5, -Tr_7-19916104703363794436192602572855262757233651585921233511204327708566528984599582067608796513614847556964923600683235803960874801547547196633123901275621828172341634280329114372437238306709825607526543371855091479474523703566741479440663807675536797962849929900204089576290912/4806288502277150141200802838468526951915129258214603774801416509601051251611778125, -b_1, -b_2, -b_3, -b_4, -b_5, -d_1, -d_2, -d_3, -d_4, -d_5, N*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, d_5, b_5, S_5, Tr_4, In_4, d_4, b_4, S_4, Tr_3, In_3, d_3, b_3, S_3, Tr_2, In_2, d_2, b_2, S_2, Tr_1, In_1, d_1, b_1, S_1, Tr_0, In_0, d_0, b_0, S_0, z_aux, w_aux, N, a, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [b, d] 252
# 89, [1 = 9, S_3*Tr_1 = 2, S_1*In_0 = 2, S_0*Tr_0 = 2, Tr_5 = 2, S_2*In_0 = 2, S_2*Tr_1 = 2, In_0*S_4 = 2, Tr_1 = 2, S_4*Tr_0 = 2, S_1*In_1 = 2, S_3*Tr_0 = 2, In_0*S_0 = 2, S_0*Tr_1 = 2, S_0*In_1 = 2, In_4*S_0 = 2, In_1*S_2 = 2, Tr_7 = 2, S_1*Tr_0 = 2, S_0*Tr_2 = 2, In_3*S_1 = 2, Tr_1*S_1 = 2, In_2*S_2 = 2, In_2*S_0 = 2, In_2*S_1 = 2, S_1*Tr_2 = 2, S_3*In_0 = 2, S_1*Tr_3 = 2, Tr_2 = 2, Tr_3*S_0 = 2, S_0*Tr_4 = 2, Tr_6 = 2, Tr_3 = 2, In_3*S_0 = 2, In_1*S_3 = 2, Tr_4 = 2, S_2*Tr_2 = 2, S_2*Tr_0 = 2, In_0*g = 1, g*In_2 = 1, In_3*S_2 = 1, g*In_3 = 1, N*S_2 = 1, In_4*N*a = 1, In_2*N*g = 1, S_1*Tr_4 = 1, g*In_5 = 1, N*In_2 = 1, In_2*S_3 = 1, S_3*Tr_2 = 1, In_4*N*g = 1, In_1*S_4 = 1, In_2*N*a = 1, Tr_0 = 1, In_1*N*a = 1, N*In_5 = 1, N*S_5 = 1, g*In_6 = 1, In_0*S_5 = 1, N*In_6 = 1, S_2*Tr_3 = 1, g*In_4 = 1, z_aux*N = 1, nu*Tr_1 = 1, In_5*N*g = 1, N*S_3 = 1, In_1*N*g = 1, In_0*N*a = 1, In_3*N*g = 1, In_1*N = 1, In_5*S_0 = 1, S_4*Tr_1 = 1, nu*Tr_2 = 1, Tr_0*nu = 1, In_5*N*a = 1, In_4*S_1 = 1, g*In_1 = 1, N*S_4 = 1, N*S_1 = 1, In_0*N*g = 1, In_3*N*a = 1, nu*Tr_6 = 1, S_5*Tr_0 = 1, nu*Tr_5 = 1, S_0*Tr_5 = 1, nu*Tr_3 = 1, N*In_3 = 1, N*In_4 = 1, nu*Tr_4 = 1]
# 134, -4.367482244
# 38
# [b = 72, d = 0]
# [b = [4, 3, 4, 4, 3, 3, 4, 3, 4, 4, 4, 3, 3, 3, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3], d = []]