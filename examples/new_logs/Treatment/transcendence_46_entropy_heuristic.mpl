infolevel[Groebner]:=10;
Et_hat := [27134232940557-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 4317375725518-b_0, b_1, -Tr_1-461728766594729075062015543, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b_0*d+In_0*N*a+In_0*N*g-In_0*S_0*b_0+In_1*N, -Tr_2+1124402428448392662985209256135281467786362940339814601108750587/65653416499, -In_2*g+Tr_2*nu+Tr_3, ((-Tr_1*d-In_1)*b_0-b_1*(Tr_0*d+In_0))*S_0-S_1*(Tr_0*d+In_0)*b_0+N*((a+g)*In_1+In_2), S_0*Tr_0*b_0*d+In_0*S_0*b_0+N*S_1, -Tr_3-91281336486453980619282010948856812139634011339357131762208898021485631171075128552296411849699479543/107759277449779135425025, -In_3*g+Tr_3*nu+Tr_4, ((-Tr_2*d-In_2)*b_0+(-Tr_0*b_2-2*Tr_1*b_1)*d-In_0*b_2-2*In_1*b_1)*S_0+((-2*S_1*Tr_1-S_2*Tr_0)*d-S_2*In_0-2*S_1*In_1)*b_0-2*S_1*Tr_0*b_1*d-2*In_0*S_1*b_1+((a+g)*In_2+In_3)*N, ((Tr_1*d+In_1)*b_0+b_1*(Tr_0*d+In_0))*S_0+S_1*(Tr_0*d+In_0)*b_0+N*S_2, b_2, -Tr_4+7410409458333446358151897451743635980999952854132832962849384098858957507565155938662255980740260934807747731664121641201302085221840376419/176869118101041203340482292812186875, -In_4*g+Tr_4*nu+Tr_5, ((-Tr_0*b_3-3*Tr_1*b_2-3*Tr_2*b_1-Tr_3*b_0)*S_0+(-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*b_0+(-3*Tr_0*b_2-6*Tr_1*b_1)*S_1-3*Tr_0*S_2*b_1)*d+(-In_0*b_3-3*In_1*b_2-3*In_2*b_1-In_3*b_0)*S_0+(-In_0*S_3-3*In_1*S_2-3*In_2*S_1)*b_0+(-3*In_0*b_2-6*In_1*b_1)*S_1-3*In_0*S_2*b_1+((a+g)*In_3+In_4)*N, ((Tr_2*d+In_2)*b_0+(Tr_0*b_2+2*Tr_1*b_1)*d+In_0*b_2+2*In_1*b_1)*S_0+((2*S_1*Tr_1+S_2*Tr_0)*d+S_2*In_0+2*S_1*In_1)*b_0+2*S_1*Tr_0*b_1*d+2*In_0*S_1*b_1+N*S_3, b_3, -Tr_5-601592510078581474899683058658452919601517231616634038484885575603841622904774974105671828980438275362242165956523758220557469147432435105437092318688194202784029691967702090927/290301546912461952211820851938324472184906265625, -In_5*g+Tr_5*nu+Tr_6, ((-Tr_0*b_4-4*Tr_1*b_3-6*Tr_2*b_2-4*Tr_3*b_1-Tr_4*b_0)*S_0+(-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*b_0+(-4*Tr_0*b_3-12*Tr_1*b_2-12*Tr_2*b_1)*S_1+(-12*S_2*Tr_1-4*S_3*Tr_0)*b_1-6*b_2*Tr_0*S_2)*d+(-In_0*b_4-4*In_1*b_3-6*In_2*b_2-4*In_3*b_1-In_4*b_0)*S_0+(-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1)*b_0+(-4*In_0*b_3-12*In_1*b_2-12*In_2*b_1)*S_1+(-4*In_0*S_3-12*In_1*S_2)*b_1-6*In_0*S_2*b_2+((a+g)*In_4+In_5)*N, ((Tr_0*b_3+3*Tr_1*b_2+3*Tr_2*b_1+Tr_3*b_0)*S_0+(3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*b_0+(3*Tr_0*b_2+6*Tr_1*b_1)*S_1+3*Tr_0*S_2*b_1)*d+(In_0*b_3+3*In_1*b_2+3*In_2*b_1+In_3*b_0)*S_0+(In_0*S_3+3*In_1*S_2+3*In_2*S_1)*b_0+(3*In_0*b_2+6*In_1*b_1)*S_1+3*In_0*S_2*b_1+N*S_4, b_4, -Tr_6+48838535875552361148676479614805980882645161024945729292771215775712617819532187905602740732649181234330683992581229834918041561585757138911887127309818025727706481290338501525202729315159905989953294771432813284891/476482209243696301051332716586995700813999789958821288671875, -In_6*g+Tr_6*nu+Tr_7, ((-Tr_0*b_5-5*Tr_1*b_4-10*Tr_2*b_3-10*Tr_3*b_2-5*Tr_4*b_1-Tr_5*b_0)*S_0+(-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*b_0+(-5*Tr_0*b_4-20*Tr_1*b_3-30*Tr_2*b_2-20*Tr_3*b_1)*S_1+(-30*S_2*Tr_2-20*S_3*Tr_1-5*S_4*Tr_0)*b_1+(-10*Tr_0*b_3-30*Tr_1*b_2)*S_2-10*b_2*Tr_0*S_3)*d+(-In_0*b_5-5*In_1*b_4-10*In_2*b_3-10*In_3*b_2-5*In_4*b_1-In_5*b_0)*S_0+(-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1)*b_0+(-5*In_0*b_4-20*In_1*b_3-30*In_2*b_2-20*In_3*b_1)*S_1+(-5*In_0*S_4-20*In_1*S_3-30*In_2*S_2)*b_1+(-10*In_0*b_3-30*In_1*b_2)*S_2-10*In_0*S_3*b_2+N*((a+g)*In_5+In_6), ((Tr_0*b_4+4*Tr_1*b_3+6*Tr_2*b_2+4*Tr_3*b_1+Tr_4*b_0)*S_0+(4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*b_0+(4*Tr_0*b_3+12*Tr_1*b_2+12*Tr_2*b_1)*S_1+(12*S_2*Tr_1+4*S_3*Tr_0)*b_1+6*b_2*Tr_0*S_2)*d+(In_0*b_4+4*In_1*b_3+6*In_2*b_2+4*In_3*b_1+In_4*b_0)*S_0+(In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1)*b_0+(4*In_0*b_3+12*In_1*b_2+12*In_2*b_1)*S_1+(4*In_0*S_3+12*In_1*S_2)*b_1+6*In_0*S_2*b_2+N*S_5, b_5, -Tr_7-3964814299574819054415867112598317356628809214822279067113284491551825069122562323061108604229099630546340869723588047649351971815317469763880824745690350984070421088294873095006430842137129677001265451409461894942395862888237083315418784613722013572103/782067123446001526079720960527778862816598038506625303107062998056640625, -b_1, -b_2, -b_3, -b_4, -b_5, -2981843798989021443560239204002012523518812132047844472322932234741663953009593617991854/107759277449779135425025-In_2, N*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, b_5, S_5, Tr_4, In_4, b_4, S_4, Tr_3, In_3, b_3, S_3, Tr_2, In_2, b_2, S_2, Tr_1, In_1, b_1, S_1, Tr_0, In_0, b_0, S_0, z_aux, w_aux, N, a, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [In_2, b] 262
# 130, [1 = 5, S_0*Tr_2*d = 1, -S_2*Tr_0*d = 1, -5*S_1*Tr_4*d = 1, -S_4*Tr_0*d = 1, In_3*N*a = 1, -Tr_6 = 1, z_aux*N = 1, -In_5*S_0 = 1, Tr_5 = 1, -4*S_3*Tr_1*d = 1, In_4*S_0 = 1, S_3*Tr_1*d = 1, -Tr_0*S_0*d = 1, -2*S_1*Tr_1*d = 1, N*In_5 = 1, -g*In_1 = 1, Tr_6 = 1, -In_0*S_4 = 1, -Tr_1 = 1, -S_2*In_0 = 1, nu*Tr_5 = 1, In_5*N*a = 1, S_2*Tr_2*d = 1, -S_0*Tr_5*d = 1, S_2*Tr_1*d = 1, In_3*S_0 = 1, -g = 1, Tr_4 = 1, S_1*Tr_1*d = 1, Tr_2 = 1, -3*In_1*S_2 = 1, N*S_4 = 1, In_1*S_2 = 1, -S_1*In_0 = 1, S_0*Tr_4*d = 1, -4*S_1*Tr_3*d = 1, In_4*N*a = 1, N*S_3 = 1, S_1*Tr_2*d = 1, -1315682618265560940163313046828858011333812427024822385846908360034690420890578284614103536152227986515/4429162478101181369991168 = 1, S_4*Tr_0*d = 1, -5*S_4*Tr_1*d = 1, In_0*S_4 = 1, -2*S_1*In_1 = 1, Tr_3 = 1, -137131846438897934395989488259513356919440535918040442735158691448033597722127124514030012395948403003318027244249234312097566692183300859744590976151210030006726634881317777196699/941639052356131101066492991853851086349459104202752 = 1, -S_3*In_0 = 1, N*In_1 = 1, -10*S_3 = 1, -10*S_2*Tr_3*d = 1, S_1*Tr_0*d = 1, -g*In_5 = 1, nu*Tr_3 = 1, -S_1*Tr_0*d = 1, nu*Tr_4 = 1, S_2 = 1, Tr_7 = 1, -S_3*Tr_0*d = 1, -Tr_1*d*S_0 = 1, -In_4*S_0 = 1, S_1*Tr_3*d = 1, -6*S_2*Tr_2*d = 1, -S_0 = 1, -10*S_3*Tr_2*d = 1, N*S_5 = 1, S_2*Tr_0*d = 1, N*S_2 = 1, -6*S_2 = 1, nu*Tr_1 = 1, In_1*S_3 = 1, In_4*N*g = 1, N*In_4 = 1, N*a = 1, -3*S_1 = 1, In_1*N*a = 1, -S_5*Tr_0*d = 1, -In_0*S_0 = 1, N*S_1 = 1, S_0*In_1 = 1, S_0*Tr_3*d = 1, -S_0*Tr_4*d = 1, S_0 = 1, -Tr_5 = 1, In_0*N*a = 1, N*g = 1, -S_0*In_1 = 1, -In_3*S_0 = 1, -5*In_1*S_4 = 1, nu*Tr_2 = 1, -In_0*g = 1, Tr_1 = 1, -14293069656421177129668895666289726602298757009023824704035425691435528055925800284357744673905162593275692842469811610820138937471441223381929361681496969335598060646325494757719137528504411136647932366518911763806099843679296216581482549850749050994510691/200192273213305425381461522839878900652085495352916541801617189551452222128128 = 1, N*In_6 = 1, Tr_1*d*S_0 = 1, -Tr_7 = 1, S_3*In_0 = 1, -Tr_2 = 1, -In_0*S_5 = 1, -10*In_3*S_2 = 1, -4*In_1*S_3 = 1, -g*In_3 = 1, In_1*N*g = 1, -S_0*Tr_3*d = 1, In_0*S_0 = 1, In_5*N*g = 1, Tr_0*S_0*d = 1, Tr_0*nu = 1, S_3*Tr_0*d = 1, N*In_3 = 1, In_3*N*g = 1, -S_0*Tr_2*d = 1, -Tr_3 = 1, S_1*In_1 = 1, In_3*S_1 = 1, -3*S_1*Tr_2*d = 1, -g*In_4 = 1, In_0*N*g = 1, S_1 = 1, S_1*In_0 = 1, N = 1, -4*In_3*S_1 = 1, -3*S_2*Tr_1*d = 1, -Tr_4 = 1, -g*In_6 = 1, -Tr_0 = 1, nu*Tr_6 = 1, -5*In_4*S_1 = 1, -1 = 1, S_2*In_0 = 1]
# 134, -4.837786147
# 1
# [In_2 = 11, b = 65]
# [In_2 = [2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3], b = [4, 3, 4, 4, 3, 3, 4, 3, 4, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3]]