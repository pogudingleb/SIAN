infolevel[Groebner]:=10;
Et_hat := [5271122508016-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 30624076928973-N_0, N_1, -Tr_1+145700644832716937752755275, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b*d+In_0*N_0*a+In_0*N_0*g-In_0*S_0*b+In_1*N_0, -Tr_2+6015910949285093070297342480571529779069756656587559916338161000/928002331181, -In_2*g+Tr_2*nu+Tr_3, (-S_0*b+(a+g)*N_0+N_1)*In_1+(-S_1*In_0-d*(S_0*Tr_1+S_1*Tr_0))*b+N_1*(a+g)*In_0+N_0*In_2, S_0*Tr_0*b*d+In_0*S_0*b+N_0*S_1, -Tr_3-16054504589894788154029222502107197130191645185212015111617772934735334854387631872469947792132042004500/9473071593451074453402371, -In_3*g+Tr_3*nu+Tr_4, (-2*S_1*In_1-In_2*S_0-S_2*In_0-(S_0*Tr_2+2*S_1*Tr_1+S_2*Tr_0)*d)*b+((2*a+2*g)*N_1+N_2)*In_1+((a+g)*N_0+2*N_1)*In_2+N_2*(a+g)*In_0+N_0*In_3, N_2, ((Tr_1*d+In_1)*S_0+S_1*(Tr_0*d+In_0))*b+N_0*S_2+S_1*N_1, -Tr_4+42844237522676956529568988491844525172134891651307991359101544701248224704835639209817967682929999389622083891506350232564936655862902810920250/96701357743838181241892552094918931661, -In_4*g+Tr_4*nu+Tr_5, ((-S_0*Tr_3-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*d-S_3*In_0-3*In_1*S_2-3*In_2*S_1-In_3*S_0)*b+(3*N_2*a+3*N_2*g+N_3)*In_1+(3*N_1*a+3*N_1*g+3*N_2)*In_2+(N_0*a+N_0*g+3*N_1)*In_3+In_0*N_3*a+In_0*N_3*g+N_0*In_4, N_3, ((2*Tr_1*d+2*In_1)*S_1+(Tr_0*d+In_0)*S_2+S_0*(Tr_2*d+In_2))*b+N_0*S_3+2*S_2*N_1+S_1*N_2, -Tr_5-114337298832298261986960384008673342152776805950960073761743188222261671699051923155785903839939761241970237454006410055597594228715487705491188386868143765866333976785395512090669250/987129939561146467140026214004468791432798312638051, -In_5*g+Tr_5*nu+Tr_6, ((-S_0*Tr_4-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*d-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1-In_4*S_0)*b+(In_0*N_4+4*In_1*N_3+6*In_2*N_2+4*In_3*N_1+In_4*N_0)*a+(In_0*N_4+4*In_1*N_3+6*In_2*N_2+4*In_3*N_1+In_4*N_0)*g+In_1*N_4+4*In_2*N_3+6*In_3*N_2+4*In_4*N_1+N_0*In_5, N_4, ((S_0*Tr_3+3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*d+3*In_1*S_2+3*In_2*S_1+In_3*S_0+S_3*In_0)*b+3*S_3*N_1+N_0*S_4+3*S_2*N_2+S_1*N_3, -Tr_6+305128966231336329806130638724624646675696407262463533919951687501109104886160799729526586577822930723646239458064686632021135773748361103926297671844431642867180614846321190074745223088035183130390845870459082634889027250/10076647736004339133588194146045562903613183667533460400728050541, -In_6*g+Tr_6*nu+Tr_7, ((-S_0*Tr_5-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*d-5*In_1*S_4-In_0*S_5-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1-In_5*S_0)*b+(In_0*N_5+5*In_1*N_4+10*In_2*N_3+10*In_3*N_2+5*In_4*N_1+In_5*N_0)*a+(In_0*N_5+5*In_1*N_4+10*In_2*N_3+10*In_3*N_2+5*In_4*N_1+In_5*N_0)*g+In_1*N_5+5*In_2*N_4+10*In_3*N_3+10*In_4*N_2+5*In_5*N_1+N_0*In_6, N_5, ((S_0*Tr_4+4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*d+In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1+In_4*S_0)*b+4*S_2*N_3+6*S_3*N_2+4*S_4*N_1+N_0*S_5+S_1*N_4, -Tr_7-74026337372154934158941679138464474904466851166620760594593679301602404069864205439657138104345079233519398399728140836887671502489290228824897435583154425030785191023323775827551986673042998741544887390194397470229138105208014045241551223415666431253731230750/9351152589501772582301149944790299947194409651356166516195381331665588218921, -N_1, -N_2, -N_3, -N_4, -N_5, -12714534866681388947045994182193746826582461934056625/10208025642991-S_1, N_0*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, S_5, N_5, Tr_4, In_4, S_4, N_4, Tr_3, In_3, S_3, N_3, Tr_2, In_2, S_2, N_2, Tr_1, In_1, S_1, N_1, Tr_0, In_0, S_0, N_0, z_aux, w_aux, a, b, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [N, S_1] 123
# 129, [1 = 6, Tr_3*b*d = 1, -10*In_2*S_3*b = 1, -2*Tr_1*d*b = 1, In_2*S_0*b = 1, S_0*Tr_4*b*d = 1, g*In_2 = 1, g*In_1 = 1, -Tr_6 = 1, S_4*Tr_0*b*d = 1, Tr_5 = 1, -5*In_1*S_4*b = 1, S_0*Tr_2*b*d = 1, In_0*b = 1, -2*b*In_1 = 1, -In_0*S_4*b = 1, -3*Tr_2*b*d = 1, -4*In_1*S_3*b = 1, -g*In_1 = 1, Tr_6 = 1, In_0*a = 1, -Tr_1 = 1, In_6 = 1, -5*Tr_4*b*d = 1, nu*Tr_5 = 1, g*In_3 = 1, -S_2*Tr_0*b*d = 1, In_1*S_0*b = 1, S_4 = 1, -In_0*S_5*b = 1, In_2*b = 1, -b*S_0*In_0 = 1, Tr_4 = 1, Tr_2 = 1, -3*In_2*b = 1, -6*In_2*S_2*b = 1, -S_5*Tr_0*b*d = 1, In_4*a = 1, S_3 = 1, In_3*a = 1, -In_3*S_0*b = 1, b*In_1 = 1, -1315682618265560940163313046828858011333812427024822385846908360034690420890578284614103536152227986515/4429162478101181369991168 = 1, -In_0*S_2*b = 1, Tr_3 = 1, S_0*Tr_1*b*d = 1, -137131846438897934395989488259513356919440535918040442735158691448033597722127124514030012395948403003318027244249234312097566692183300859744590976151210030006726634881317777196699/941639052356131101066492991853851086349459104202752 = 1, Tr_0*b*d = 1, S_0*Tr_3*b*d = 1, Tr_1*d*b = 1, -In_2*S_0*b = 1, -5*S_4*Tr_1*b*d = 1, -g*In_5 = 1, nu*Tr_3 = 1, nu*Tr_4 = 1, -g*In_2 = 1, In_2 = 1, S_5 = 1, -In_4*S_0*b = 1, S_2 = 1, Tr_7 = 1, In_3 = 1, -10*S_3*Tr_2*b*d = 1, Tr_2*b*d = 1, g*In_4 = 1, -4*In_3*b = 1, -S_0*Tr_4*b*d = 1, -3*In_1*S_2*b = 1, In_4 = 1, In_1*S_3*b = 1, g*In_5 = 1, In_1*S_2*b = 1, nu*Tr_1 = 1, d*b*S_0*Tr_0 = 1, In_3*S_0*b = 1, S_2*Tr_0*b*d = 1, S_3*Tr_0*b*d = 1, In_0*S_2*b = 1, -S_0*Tr_2*b*d = 1, In_0*g = 1, -Tr_5 = 1, -4*Tr_3*b*d = 1, nu*Tr_2 = 1, -In_5*S_0*b = 1, -In_0*g = 1, Tr_1 = 1, S_2*Tr_2*b*d = 1, -10*S_2*Tr_3*b*d = 1, -In_1*S_0*b = 1, -14293069656421177129668895666289726602298757009023824704035425691435528055925800284357744673905162593275692842469811610820138937471441223381929361681496969335598060646325494757719137528504411136647932366518911763806099843679296216581482549850749050994510691/200192273213305425381461522839878900652085495352916541801617189551452222128128 = 1, In_0*S_4*b = 1, -Tr_7 = 1, -Tr_2 = 1, -g*In_3 = 1, In_1 = 1, S_2*Tr_1*b*d = 1, In_2*S_2*b = 1, b*S_0*In_0 = 1, In_5*a = 1, Tr_0*nu = 1, In_0*S_3*b = 1, -d*b*S_0*Tr_0 = 1, In_4*S_0*b = 1, In_5 = 1, -Tr_3 = 1, In_1*a = 1, -3*S_2*Tr_1*b*d = 1, In_2*a = 1, -g*In_4 = 1, -In_0*b = 1, In_3*b = 1, -6*S_2*Tr_2*b*d = 1, -In_0*S_3*b = 1, -S_0*Tr_1*b*d = 1, S_3*Tr_1*b*d = 1, -S_4*Tr_0*b*d = 1, -S_0*Tr_5*b*d = 1, z_aux = 1, -4*S_3*Tr_1*b*d = 1, -Tr_4 = 1, -g*In_6 = 1, -Tr_0*b*d = 1, -Tr_0 = 1, -5*In_4*b = 1, -S_3*Tr_0*b*d = 1, nu*Tr_6 = 1, -S_0*Tr_3*b*d = 1, -10*In_3*S_2*b = 1, -1 = 1]
# 134, -4.817611764
# 1
# [N = 24, S_1 = 18]
# [N = [3, 3, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 2], S_1 = [4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3]]