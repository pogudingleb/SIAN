infolevel[Groebner]:=10;
Et_hat := [96530368010292-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 226492464962683-N_0, N_1, 81895059920095-b_0, b_1, -Tr_1-6513918573690477180499884316, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b_0*d+In_0*N_0*a+In_0*N_0*g-In_0*S_0*b_0+In_1*N_0, -Tr_2+2085942163471587288548477779440735702898085248535100182621616836982844/226492464962683, -In_2*g+Tr_2*nu+Tr_3, (-b_1*S_0-S_1*b_0+N_1*(a+g))*In_0+(-b_0*S_0+(a+g)*N_0+N_1)*In_1-d*(Tr_0*b_1+Tr_1*b_0)*S_0-S_1*Tr_0*b_0*d+N_0*In_2, S_0*Tr_0*b_0*d+In_0*S_0*b_0+N_0*S_1, -Tr_3-2793882100023223344248996749990827055869983210904430308848106887972123155578108631155634609385412349823898759996/51298836684872186368582558489, -In_3*g+Tr_3*nu+Tr_4, ((-Tr_2*d-In_2)*b_0+(-Tr_0*b_2-2*Tr_1*b_1)*d-2*b_1*In_1-b_2*In_0)*S_0+((-2*S_1*Tr_1-S_2*Tr_0)*d-S_2*In_0-2*S_1*In_1)*b_0-2*S_1*Tr_0*b_1*d+(-2*b_1*S_1+N_2*(a+g))*In_0+((2*a+2*g)*N_1+N_2)*In_1+((a+g)*N_0+2*N_1)*In_2+N_0*In_3, N_2, ((Tr_1*d+In_1)*b_0+b_1*(Tr_0*d+In_0))*S_0+S_1*(Tr_0*d+In_0)*b_0+N_0*S_2+S_1*N_1, b_2, -Tr_4+3742087065270661720764858737657090558892990802640832473205741247960148206526883169575138924097546015916717221968516057155325963992597679132854931546007964/11618799970474811011990286849463890049865987, -In_4*g+Tr_4*nu+Tr_5, ((-Tr_0*b_3-3*Tr_1*b_2-3*Tr_2*b_1-Tr_3*b_0)*S_0+(-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*b_0+(-3*Tr_0*b_2-6*Tr_1*b_1)*S_1-3*b_1*Tr_0*S_2)*d+(-In_0*b_3-3*In_1*b_2-3*In_2*b_1-In_3*b_0)*S_0+(-In_0*S_3-3*In_1*S_2-3*In_2*S_1)*b_0+(-3*b_2*S_1-3*S_2*b_1+N_3*(a+g))*In_0+(3*N_2*a+3*N_2*g-6*S_1*b_1+N_3)*In_1+(3*N_1*a+3*N_1*g+3*N_2)*In_2+(N_0*a+N_0*g+3*N_1)*In_3+N_0*In_4, N_3, ((Tr_2*d+In_2)*b_0+(Tr_0*b_2+2*Tr_1*b_1)*d+2*b_1*In_1+b_2*In_0)*S_0+((2*S_1*Tr_1+S_2*Tr_0)*d+S_2*In_0+2*S_1*In_1)*b_0+2*S_1*Tr_0*b_1*d+(2*In_0*b_1+N_2)*S_1+N_0*S_3+2*S_2*N_1, b_3, -Tr_5-5012099688798397087489445066172915754226338337749786719614505189031758532440767314563149677641027429197556931185233775251955376610775444310198660439197016801176623142334916158523515828411034522076/2631570645221188408016616102057718851521676330297105963121, -In_5*g+Tr_5*nu+Tr_6, ((-Tr_0*b_4-4*Tr_1*b_3-6*Tr_2*b_2-4*Tr_3*b_1-Tr_4*b_0)*S_0+(-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*b_0+(-4*Tr_0*b_3-12*Tr_1*b_2-12*Tr_2*b_1)*S_1+(-12*S_2*Tr_1-4*S_3*Tr_0)*b_1-6*b_2*Tr_0*S_2)*d+(-In_0*b_4-4*In_1*b_3-6*In_2*b_2-4*In_3*b_1-In_4*b_0)*S_0+(-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1)*b_0+(-4*In_0*b_3-12*In_1*b_2-12*In_2*b_1)*S_1+(-4*In_0*S_3-12*In_1*S_2)*b_1+(-6*S_2*b_2+N_4*(a+g))*In_0+(4*N_3*a+4*N_3*g+N_4)*In_1+(6*N_2*a+6*N_2*g+4*N_3)*In_2+(4*N_1*a+4*N_1*g+6*N_2)*In_3+In_4*N_0*a+In_4*N_0*g+4*In_4*N_1+N_0*In_5, N_4, ((Tr_0*b_3+3*Tr_1*b_2+3*Tr_2*b_1+Tr_3*b_0)*S_0+(3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*b_0+(3*Tr_0*b_2+6*Tr_1*b_1)*S_1+3*b_1*Tr_0*S_2)*d+(In_0*b_3+3*In_1*b_2+3*In_2*b_1+In_3*b_0)*S_0+(In_0*S_3+3*In_1*S_2+3*In_2*S_1)*b_0+(3*In_0*b_2+6*In_1*b_1+N_3)*S_1+3*In_0*S_2*b_1+N_0*S_4+3*S_3*N_1+3*S_2*N_2, b_4, -Tr_6+6713137041517531349685925416297051434189422598973002065969307647963030051175297069212628900497652239710842371299135910986776122116565313212412321212501907724425492337454517070985720735189284007807187832444187212353172225687463239885730684/596030922159585110993390054091788251877625578598911316335777280039213643, -In_6*g+Tr_6*nu+Tr_7, ((-Tr_0*b_5-5*Tr_1*b_4-10*Tr_2*b_3-10*Tr_3*b_2-5*Tr_4*b_1-Tr_5*b_0)*S_0+(-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*b_0+(-5*Tr_0*b_4-20*Tr_1*b_3-30*Tr_2*b_2-20*Tr_3*b_1)*S_1+(-30*S_2*Tr_2-20*S_3*Tr_1-5*S_4*Tr_0)*b_1+(-10*Tr_0*b_3-30*Tr_1*b_2)*S_2-10*S_3*Tr_0*b_2)*d+(-In_0*b_5-5*In_1*b_4-10*In_2*b_3-10*In_3*b_2-5*In_4*b_1-In_5*b_0)*S_0+(-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1)*b_0+(-5*In_0*b_4-20*In_1*b_3-30*In_2*b_2-20*In_3*b_1)*S_1+(-5*In_0*S_4-20*In_1*S_3-30*In_2*S_2)*b_1+(-10*S_2*b_3-10*S_3*b_2+N_5*(a+g))*In_0+(5*N_4*a+5*N_4*g-30*S_2*b_2+N_5)*In_1+(10*N_3*a+10*N_3*g+5*N_4)*In_2+(10*N_2*a+10*N_2*g+10*N_3)*In_3+(5*In_4*N_1+In_5*N_0)*a+(5*In_4*N_1+In_5*N_0)*g+10*In_4*N_2+5*In_5*N_1+N_0*In_6, N_5, ((Tr_0*b_4+4*Tr_1*b_3+6*Tr_2*b_2+4*Tr_3*b_1+Tr_4*b_0)*S_0+(4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*b_0+(4*Tr_0*b_3+12*Tr_1*b_2+12*Tr_2*b_1)*S_1+(12*S_2*Tr_1+4*S_3*Tr_0)*b_1+6*b_2*Tr_0*S_2)*d+(In_0*b_4+4*In_1*b_3+6*In_2*b_2+4*In_3*b_1+In_4*b_0)*S_0+(In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1)*b_0+(4*In_0*b_3+12*In_1*b_2+12*In_2*b_1+N_4)*S_1+(4*In_0*S_3+12*In_1*S_2)*b_1+(6*In_0*b_2+4*N_3)*S_2+6*S_3*N_2+4*S_4*N_1+N_0*S_5, b_5, -Tr_7-8991482958511426819521762137615964564303020029507552092314773469037122206665858932157991281306979954129399023034325966045193129179021903504272006861798158787392136212028980455609609969910554703917584759282413180343968602411210983968444245665471720006507045061337446347191073442556/134996512753905469243962274470792120777261033448601316636709390255374593905476659484169, -N_1, -N_2, -N_3, -N_4, -N_5, -b_1, -b_2, -b_3, -b_4, -b_5, N_0*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, b_5, S_5, N_5, Tr_4, In_4, b_4, S_4, N_4, Tr_3, In_3, b_3, S_3, N_3, Tr_2, In_2, b_2, S_2, N_2, Tr_1, In_1, b_1, S_1, N_1, Tr_0, In_0, b_0, S_0, N_0, z_aux, w_aux, a, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [N, b] 312
# 131, [1 = 4, -Tr_4 = 1, -10*In_2*S_3 = 1, In_0*g = 1, -S_0*Tr_3*d = 1, S_1*In_1 = 1, g*In_3 = 1, g*In_1 = 1, S_3*Tr_0*d = 1, In_2 = 1, -179309634593493733658162212844461858518738205059942771396968370122660798111006850172798173061927452435076970925532086626741233547613316281809745770255690419192793534461037590319412407502817305399172689471404860699487289686616281977054599048429017935279549578831136/246035506446721274012617103349915219301813950982064232570241051996976213527 = 1, -195729598087814495425614884200452077342696811229758649392798447939034136730600123569559915378863964696973980549415713899267784372098356115932139834648691705675571488195208248627151963570/274850520894250566433945179777030662974823983471063 = 1, -In_3*S_0 = 1, -3*S_1*Tr_2*d = 1, -5*S_4*Tr_1*d = 1, S_1*In_0 = 1, -Tr_3 = 1, -S_5*Tr_0*d = 1, g*In_4 = 1, -d*Tr_0*S_0 = 1, In_1*a = 1, -5*In_1*S_4 = 1, S_4 = 1, S_1*Tr_3*d = 1, S_3 = 1, Tr_5 = 1, S_2*Tr_1*d = 1, -Tr_6 = 1, -In_0*S_5 = 1, In_3*S_0 = 1, S_2*Tr_0*d = 1, In_2*S_0 = 1, In_2*S_2 = 1, In_4*a = 1, -S_1*In_0 = 1, -Tr_7 = 1, -g*In_1 = 1, Tr_7 = 1, -S_3*Tr_0*d = 1, In_0*S_4 = 1, Tr_3 = 1, Tr_4 = 1, -5*S_1*Tr_4*d = 1, In_1 = 1, nu*Tr_5 = 1, Tr_0*nu = 1, nu*Tr_3 = 1, -S_0*Tr_4*d = 1, -S_0*In_1 = 1, In_3 = 1, In_0*a = 1, -In_5*g = 1, -5*In_4*S_1 = 1, nu*Tr_4 = 1, -4360269199474597815095572546544347080173935020453374789825453251201657350086524463723247049405869120527092/6266128001226122488850503 = 1, -In_5*S_0 = 1, Tr_6 = 1, -10*S_2*Tr_3*d = 1, S_3*In_0 = 1, -6*In_2*S_2 = 1, -S_1*Tr_0*d = 1, In_1*S_3 = 1, nu*Tr_2 = 1, S_4*Tr_0*d = 1, In_6 = 1, S_1 = 1, -g*In_4 = 1, S_0*Tr_4*d = 1, In_4 = 1, -Tr_1*d*S_0 = 1, -10*In_3*S_2 = 1, -Tr_2 = 1, -10*S_3*Tr_2*d = 1, -4*S_1*Tr_3*d = 1, S_5 = 1, Tr_1*d*S_0 = 1, -Tr_0 = 1, -4*In_3*S_1 = 1, -3*S_2*Tr_1*d = 1, S_0*In_1 = 1, -S_0*Tr_2*d = 1, Tr_2 = 1, S_3*Tr_1*d = 1, -3*In_1*S_2 = 1, In_5 = 1, -4*In_1*S_3 = 1, -In_4*S_0 = 1, nu*Tr_6 = 1, -S_3*In_0 = 1, -2*S_1*In_1 = 1, In_2*a = 1, -3*In_2*S_1 = 1, -S_2*In_0 = 1, -S_0*Tr_5*d = 1, -g*In_2 = 1, -In_0*g = 1, -2*S_1*Tr_1*d = 1, -In_0*S_4 = 1, -S_2*Tr_0*d = 1, S_0*Tr_3*d = 1, Tr_1 = 1, S_0*Tr_2*d = 1, S_1*Tr_0*d = 1, -g*In_3 = 1, In_3*a = 1, -5429305681874946826554610 = 1, z_aux = 1, In_2*S_1 = 1, nu*Tr_1 = 1, S_0*In_0 = 1, S_2*Tr_2*d = 1, -S_4*Tr_0*d = 1, d*Tr_0*S_0 = 1, In_4*S_0 = 1, -4*S_3*Tr_1*d = 1, S_2 = 1, S_1*Tr_2*d = 1, In_5*a = 1, S_1*Tr_1*d = 1, -In_6*g = 1, In_5*g = 1, In_3*S_1 = 1, g*In_2 = 1, -S_0*In_0 = 1, S_2*In_0 = 1, -Tr_5 = 1, -Tr_1 = 1, In_1*S_2 = 1, -6*S_2*Tr_2*d = 1, -1 = 1, -In_2*S_0 = 1]
# 134, -4.856457879
# 1
# [N = 24, b = 72]
# [N = [3, 3, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 2], b = [4, 3, 4, 4, 3, 3, 4, 3, 4, 4, 4, 3, 3, 3, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3]]