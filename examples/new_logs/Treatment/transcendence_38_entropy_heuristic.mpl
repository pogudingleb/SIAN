infolevel[Groebner]:=10;
Et_hat := [9424041583660-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 23476751511557-N_0, N_1, -Tr_1+47031942743589921386271685, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b*d+In_0*N_0*a+In_0*N_0*g-In_0*S_0*b+In_1*N_0, -Tr_2+163180563551415311644299894879835362800544400695620418945848691445/23476751511557, -In_2*g+Tr_2*nu+Tr_3, (-S_0*b+(a+g)*N_0+N_1)*In_1+(-S_1*In_0-d*(S_0*Tr_1+S_1*Tr_0))*b+N_1*(a+g)*In_0+N_0*In_2, S_0*Tr_0*b*d+In_0*S_0*b+N_0*S_1, -Tr_3-161466091396632995082915710862956059029941004304268562459160414883840713881358958450179566441443413819235/551157861535393884304564249, -In_3*g+Tr_3*nu+Tr_4, (-2*S_1*In_1-In_2*S_0-S_2*In_0-(S_0*Tr_2+2*S_1*Tr_1+S_2*Tr_0)*d)*b+((2*a+2*g)*N_1+N_2)*In_1+((a+g)*N_0+2*N_1)*In_2+N_2*(a+g)*In_0+N_0*In_3, N_2, ((Tr_1*d+In_1)*S_0+S_1*(Tr_0*d+In_0))*b+N_0*S_2+S_1*N_1, -Tr_4+159769632506796277701307261322362106368505226275850607183095435357413967816702998822377625617438398905008170654549978081150624843137154640220005/12939396158907582082202552310464972525693, -In_4*g+Tr_4*nu+Tr_5, ((-S_0*Tr_3-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*d-S_3*In_0-3*In_1*S_2-3*In_2*S_1-In_3*S_0)*b+(3*N_2*a+3*N_2*g+N_3)*In_1+(3*N_1*a+3*N_1*g+3*N_2)*In_2+(N_0*a+N_0*g+3*N_1)*In_3+In_0*N_3*a+In_0*N_3*g+N_0*In_4, N_3, ((2*Tr_1*d+2*In_1)*S_1+(Tr_0*d+In_0)*S_2+S_0*(Tr_2*d+In_2))*b+N_0*S_3+2*S_2*N_1+S_1*N_2, -Tr_5-158090997624206732930356257839250696437131669286011746439257446349342628181471049196579596811602565841949860921186928623807079128213527694868454514658817641575009518564706446744205315/303774988332268417418216819382551906492065613768934001, -In_5*g+Tr_5*nu+Tr_6, ((-S_0*Tr_4-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*d-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1-In_4*S_0)*b+(In_0*N_4+4*In_1*N_3+6*In_2*N_2+4*In_3*N_1+In_4*N_0)*a+(In_0*N_4+4*In_1*N_3+6*In_2*N_2+4*In_3*N_1+In_4*N_0)*g+In_1*N_4+4*In_2*N_3+6*In_3*N_2+4*In_4*N_1+N_0*In_5, N_4, ((S_0*Tr_3+3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*d+3*In_1*S_2+3*In_2*S_1+In_3*S_0+S_3*In_0)*b+3*S_3*N_1+N_0*S_4+3*S_2*N_2+S_1*N_3, -Tr_6+156429999479425526401915569662231100249677383119637960255187638802995900980963511556572590489300855043988824534748825561768613692269454891800225358287963957635110233743968974646600344622817364442769269511587881312521489445/7131649916502792607181773941866886326169613519477044259705321749557, -In_6*g+Tr_6*nu+Tr_7, ((-S_0*Tr_5-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*d-5*In_1*S_4-In_0*S_5-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1-In_5*S_0)*b+(In_0*N_5+5*In_1*N_4+10*In_2*N_3+10*In_3*N_2+5*In_4*N_1+In_5*N_0)*a+(In_0*N_5+5*In_1*N_4+10*In_2*N_3+10*In_3*N_2+5*In_4*N_1+In_5*N_0)*g+In_1*N_5+5*In_2*N_4+10*In_3*N_3+10*In_4*N_2+5*In_5*N_1+N_0*In_6, N_5, ((S_0*Tr_4+4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*d+In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1+In_4*S_0)*b+4*S_2*N_3+6*S_3*N_2+4*S_4*N_1+N_0*S_5+S_1*N_4, -Tr_7-154786452770765595568142844949067201007867567203197749849440634128221623322645869853019521939668179512565404002378021372554546768909043331056521536692980475420762592560255697227145958460733862870012917456898488557588264984853512336491670010078719071168790511235/167427972957152289179866396323584097804387568719345201484199502451309327545130249, -N_1, -N_2, -N_3, -N_4, -N_5, -6091036095333565160804203868240096190141729040572426253115519170552675328830348136658310215/551157861535393884304564249-In_2, N_0*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, S_5, N_5, Tr_4, In_4, S_4, N_4, Tr_3, In_3, S_3, N_3, Tr_2, In_2, S_2, N_2, Tr_1, In_1, S_1, N_1, Tr_0, In_0, S_0, N_0, z_aux, w_aux, a, b, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [In_2, N] 83
# 129, [1 = 6, In_3*S_1*b = 1, S_0*Tr_4*b*d = 1, -2*S_1*Tr_1*b*d = 1, g*In_1 = 1, -Tr_6 = 1, S_4*Tr_0*b*d = 1, Tr_5 = 1, In_5*a = 1, -5*In_1*S_4*b = 1, S_0*Tr_2*b*d = 1, S_0*b = 1, S_1*Tr_2*b*d = 1, -In_0*S_4*b = 1, S_2*b = 1, -4*In_1*S_3*b = 1, -g*In_1 = 1, Tr_6 = 1, -Tr_1 = 1, In_6 = 1, nu*Tr_5 = 1, g*In_3 = 1, -S_2*Tr_0*b*d = 1, In_1*S_0*b = 1, S_4 = 1, -In_0*S_5*b = 1, -g = 1, -b*S_0*In_0 = 1, Tr_4 = 1, Tr_2 = 1, In_1*a = 1, In_0*a = 1, -S_5*Tr_0*b*d = 1, S_1*Tr_3*b*d = 1, S_3 = 1, -2*In_1*S_1*b = 1, -In_3*S_0*b = 1, -1315682618265560940163313046828858011333812427024822385846908360034690420890578284614103536152227986515/4429162478101181369991168 = 1, -In_0*S_2*b = 1, -4*In_3*S_1*b = 1, In_0*S_1*b = 1, Tr_3 = 1, S_0*Tr_1*b*d = 1, -137131846438897934395989488259513356919440535918040442735158691448033597722127124514030012395948403003318027244249234312097566692183300859744590976151210030006726634881317777196699/941639052356131101066492991853851086349459104202752 = 1, S_0*Tr_3*b*d = 1, a = 1, -5*S_4*Tr_1*b*d = 1, -g*In_5 = 1, nu*Tr_3 = 1, nu*Tr_4 = 1, S_5 = 1, -In_4*S_0*b = 1, S_2 = 1, Tr_7 = 1, In_3 = 1, -10*S_3*b = 1, -10*S_3*Tr_2*b*d = 1, g*In_4 = 1, -S_0*Tr_4*b*d = 1, -3*In_1*S_2*b = 1, -S_1*Tr_0*b*d = 1, In_4 = 1, In_1*S_3*b = 1, g*In_5 = 1, In_1*S_2*b = 1, nu*Tr_1 = 1, d*b*S_0*Tr_0 = 1, In_3*S_0*b = 1, S_2*Tr_0*b*d = 1, -5*S_1*Tr_4*b*d = 1, S_3*Tr_0*b*d = 1, In_4*a = 1, In_0*S_2*b = 1, -S_0*Tr_2*b*d = 1, In_0*g = 1, -In_0*S_1*b = 1, S_1*Tr_0*b*d = 1, -Tr_5 = 1, nu*Tr_2 = 1, -In_5*S_0*b = 1, -In_0*g = 1, Tr_1 = 1, S_2*Tr_2*b*d = 1, -4*S_1*Tr_3*b*d = 1, -10*S_2*Tr_3*b*d = 1, -3*S_1*b = 1, -In_1*S_0*b = 1, -14293069656421177129668895666289726602298757009023824704035425691435528055925800284357744673905162593275692842469811610820138937471441223381929361681496969335598060646325494757719137528504411136647932366518911763806099843679296216581482549850749050994510691/200192273213305425381461522839878900652085495352916541801617189551452222128128 = 1, In_0*S_4*b = 1, -S_0*b = 1, -Tr_7 = 1, -Tr_2 = 1, -5*In_4*S_1*b = 1, S_1*Tr_1*b*d = 1, -g*In_3 = 1, In_1 = 1, S_2*Tr_1*b*d = 1, In_1*S_1*b = 1, b*S_0*In_0 = 1, Tr_0*nu = 1, In_0*S_3*b = 1, -d*b*S_0*Tr_0 = 1, In_4*S_0*b = 1, In_5 = 1, -Tr_3 = 1, -3*S_1*Tr_2*b*d = 1, S_1*b = 1, -3*S_2*Tr_1*b*d = 1, -g*In_4 = 1, -6*S_2*b = 1, -6*S_2*Tr_2*b*d = 1, -In_0*S_3*b = 1, S_1 = 1, -S_0*Tr_1*b*d = 1, S_3*Tr_1*b*d = 1, -S_4*Tr_0*b*d = 1, In_3*a = 1, -S_0*Tr_5*b*d = 1, z_aux = 1, -4*S_3*Tr_1*b*d = 1, g = 1, -Tr_4 = 1, -g*In_6 = 1, -Tr_0 = 1, -S_3*Tr_0*b*d = 1, nu*Tr_6 = 1, -S_0*Tr_3*b*d = 1, -10*In_3*S_2*b = 1, -1 = 1]
# 134, -4.817611764
# 1
# [N = 21, In_2 = 11]
# [N = [3, 3, 2, 3, 3, 2, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 2], In_2 = [2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3]]