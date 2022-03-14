infolevel[Groebner]:=10;
Et_hat := [378206028780908863623735-x5_0, -k7*x1_0+x5_1, 252370027472131042970088-b_0, b_1, -x5_1+115759692118459193769592747242710992756220150948, -k7*x1_1+x5_2, k3*x1_0-k4*x2_0+k7*x1_0+x1_1, -x5_2-37127916561443935189831374799576569155178057796274089531947408260007140, -k7*x1_2+x5_3, (k3+k7)*x1_1+x1_2-k4*x2_1, b_0*d*k5*x2_0+a*k5*x2_0-k5*x2_0*x3_0-k5*x2_0*x4_0-k3*x1_0+k4*x2_0-k6*x3_0-k6*x4_0+x2_1, -x5_3-35803718169091439970432149133262519642048593743911220111137691460315060552151039634632949297276298879550066406871410222973252514765886857684, -k7*x1_3+x5_4, (k3+k7)*x1_2+x1_3-k4*x2_2, ((b_0*d+a-x3_0-x4_0)*x2_1-x2_0*(-b_1*d+x3_1+x4_1))*k5+k4*x2_1-k6*x3_1-k6*x4_1-x1_1*k3+x2_2, -a*k5*x2_0+k5*x2_0*x3_0+k6*x3_0+x3_1, -b_0*d*k5*x2_0+k5*x2_0*x4_0+k6*x4_0+x4_1, -x5_4+120458300186953252892613115942130714879046925979706393644764556423115050018495057320248585604245561657532507382941793483111215298004799161207043321321954616551752698777989899383053749197135738459452984208553340, -k7*x1_4+x5_5, (k3+k7)*x1_3+x1_4-k4*x2_3, ((b_0*d+a-x3_0-x4_0)*x2_2+(b_2*d-x3_2-x4_2)*x2_0-2*x2_1*(-b_1*d+x3_1+x4_1))*k5-x1_2*k3+k4*x2_2-k6*x3_2-k6*x4_2+x2_3, b_2, ((-a+x3_0)*x2_1+x2_0*x3_1)*k5+k6*x3_1+x3_2, ((-b_1*d+x4_1)*x2_0-x2_1*(b_0*d-x4_0))*k5+k6*x4_1+x4_2, -x5_5-405270816159436743368066966557449557861666501089840832241373929817174952747745744333581253740020032236710456047489685466482963597624349761072830873082814519432927883843776716233298144251983865320404248061657169916653372939452168935197825353369807381274646002418295348083122660596, -k7*x1_5+x5_6, (k3+k7)*x1_4+x1_5-k4*x2_4, ((b_0*d+a-x3_0-x4_0)*x2_3+(3*b_1*x2_2+3*b_2*x2_1+b_3*x2_0)*d+(-x3_3-x4_3)*x2_0+(-3*x3_2-3*x4_2)*x2_1-3*x2_2*(x3_1+x4_1))*k5-x1_3*k3+k4*x2_3-k6*x3_3-k6*x4_3+x2_4, b_3, ((-a+x3_0)*x2_2+x2_0*x3_2+2*x2_1*x3_1)*k5+k6*x3_2+x3_3, ((-b_0*x2_2-2*b_1*x2_1-b_2*x2_0)*d+x4_0*x2_2+2*x2_1*x4_1+x4_2*x2_0)*k5+k6*x4_2+x4_3, -x5_6+1363496198897260915231104456238372734914209609992218659451729250917246546884756278897166219756628654827624073643409366815471150600909567163934711540530211899946807882229633443147908551057426509833935819107060506298764265479327076370985512974476488239076339215747452129718464147001849836376083378696557020917946439791911287093562584334966172356051036, -k7*x1_6+x5_7, (k3+k7)*x1_5+x1_6-k4*x2_5, ((b_0*d+a-x3_0-x4_0)*x2_4+(4*b_1*x2_3+6*b_2*x2_2+4*b_3*x2_1+b_4*x2_0)*d+(-x3_4-x4_4)*x2_0+(-4*x3_3-4*x4_3)*x2_1+(-6*x3_2-6*x4_2)*x2_2-4*x2_3*(x3_1+x4_1))*k5-x1_4*k3+k4*x2_4-k6*x3_4-k6*x4_4+x2_5, b_4, ((-a+x3_0)*x2_3+3*x3_2*x2_1+x3_3*x2_0+3*x2_2*x3_1)*k5+k6*x3_3+x3_4, ((-b_0*x2_3-3*b_1*x2_2-3*b_2*x2_1-b_3*x2_0)*d+3*x4_1*x2_2+x2_3*x4_0+3*x4_2*x2_1+x4_3*x2_0)*k5+k6*x4_3+x4_4, -x5_7-4587356923514289391408460013440284402666216855641898069608891262377956522086833019276040968803235574307946411955238909860610732741470150107652421888565504176958090625249748249931453955900998439173965723837276240192666565508193601784924019077224884150930984975395617655643484846787469708288025157814765403373895001991101320549132084805758681111195230394837865487941580496177175288579946664421168959856741753983163641364, -k7*x1_7+x5_8, (k3+k7)*x1_6+x1_7-k4*x2_6, ((b_0*x2_5+5*b_1*x2_4+10*b_2*x2_3+10*b_3*x2_2+5*b_4*x2_1+b_5*x2_0)*d+(a-x3_0-x4_0)*x2_5+(-x3_5-x4_5)*x2_0+(-5*x3_4-5*x4_4)*x2_1+(-10*x3_3-10*x4_3)*x2_2+(-10*x3_2-10*x4_2)*x2_3-5*x2_4*(x3_1+x4_1))*k5-x1_5*k3+k4*x2_5-k6*x3_5-k6*x4_5+x2_6, b_5, ((-a+x3_0)*x2_4+4*x3_1*x2_3+6*x3_2*x2_2+4*x3_3*x2_1+x3_4*x2_0)*k5+k6*x3_4+x3_5, ((-b_0*x2_4-4*b_1*x2_3-6*b_2*x2_2-4*b_3*x2_1-b_4*x2_0)*d+4*x4_3*x2_1+6*x4_2*x2_2+4*x4_1*x2_3+x4_0*x2_4+x4_4*x2_0)*k5+k6*x4_4+x4_5, -x5_8+15433738327055016659973388969738916336000402221177030706892401016168556688063178245812434413201031692058484609039613307486653526054501519886416723417598794503546207704407032369159100026818364180636291218145511154294315061876828561088607718743629819593856834539490352202607749020220975969533686587275075692236146452914604700541144628508715613476573933512947508777222643235344599584165049921073330320467507907869627493325493552010958255904441888241123506245792319481407743882666787428406844, -k7*x1_8+x5_9, (k3+k7)*x1_7+x1_8-k4*x2_7, ((b_0*x2_6+6*b_1*x2_5+15*b_2*x2_4+20*b_3*x2_3+15*b_4*x2_2+6*b_5*x2_1+b_6*x2_0)*d+(a-x3_0-x4_0)*x2_6+(-x3_6-x4_6)*x2_0+(-6*x3_5-6*x4_5)*x2_1+(-15*x3_4-15*x4_4)*x2_2+(-20*x3_3-20*x4_3)*x2_3+(-15*x3_2-15*x4_2)*x2_4-6*x2_5*(x3_1+x4_1))*k5-x1_6*k3+k4*x2_6-k6*x3_6-k6*x4_6+x2_7, b_6, ((-a+x3_0)*x2_5+5*x3_1*x2_4+10*x3_2*x2_3+10*x3_3*x2_2+5*x3_4*x2_1+x3_5*x2_0)*k5+k6*x3_5+x3_6, ((-b_0*x2_5-5*b_1*x2_4-10*b_2*x2_3-10*b_3*x2_2-5*b_4*x2_1-b_5*x2_0)*d+x2_0*x4_5+5*x4_4*x2_1+10*x4_3*x2_2+10*x4_2*x2_3+5*x4_1*x2_4+x4_0*x2_5)*k5+k6*x4_5+x4_6, -x5_9-51925385950898747901284140011119323150494482621501557692798097207439798200043165041125374551826753648652579404944825155664594365090399965661757151128117589362118929717250617397340818920261084732518390864677805393353503310030683885604901246745323610277010034494071966317002391119489057460005991397436908090781759451367621015012708551688115738851214003982328883234726483384652988924557578956702124316518075658978593938354306994561744219963263298941694605612054777904633325285445725803864444578877455562271968153611501358916391699275781128729696607856219595828, -k7*x1_9+x5_10, (k3+k7)*x1_8+x1_9-k4*x2_8, ((b_0*x2_7+7*b_1*x2_6+21*b_2*x2_5+35*b_3*x2_4+35*b_4*x2_3+21*b_5*x2_2+7*b_6*x2_1+b_7*x2_0)*d+(a-x3_0-x4_0)*x2_7+(-x3_7-x4_7)*x2_0+(-7*x3_6-7*x4_6)*x2_1+(-21*x3_5-21*x4_5)*x2_2+(-35*x3_4-35*x4_4)*x2_3+(-35*x3_3-35*x4_3)*x2_4+(-21*x3_2-21*x4_2)*x2_5-7*x2_6*(x3_1+x4_1))*k5-x1_7*k3+k4*x2_7-k6*x3_7-k6*x4_7+x2_8, b_7, ((-a+x3_0)*x2_6+6*x3_1*x2_5+15*x3_2*x2_4+20*x3_3*x2_3+15*x3_4*x2_2+6*x2_1*x3_5+x3_6*x2_0)*k5+k6*x3_6+x3_7, ((-b_0*x2_6-6*b_1*x2_5-15*b_2*x2_4-20*b_3*x2_3-15*b_4*x2_2-6*b_5*x2_1-b_6*x2_0)*d+x2_0*x4_6+6*x4_5*x2_1+15*x4_4*x2_2+20*x4_3*x2_3+15*x4_2*x2_4+6*x2_5*x4_1+x4_0*x2_6)*k5+k6*x4_6+x4_7, -x5_10+174698161198141567821614863740390123097058755299293633302214958589075947568150826710584847758712911432898220102684190965499909896863144645328799185875644871920263955411870047016331357145258184339643318927183921882059382032826627660175298539313169136989200774190435771317912871061816692666789252462597156135730363358588521948163378487813215083429465139751132505106389778948111416270150580221694474726451844874870078009412602945125267512124736265074956009239112909947756322528966567584808750118073395629392467494263126117670747684148789635111995961764513118837124432433124040662388331796485716363175785500229173619121148913235228, -b_1, -b_2, -b_3, -b_4, -b_5, -b_6, -b_7, -56954981063989016774129892726588663103125717369985392421939752108188268892170036701315126272871314571493045024892266686428492625376031169610965895245151393950616723450855954682882136485645940666312359362132243564848286185411445283650884861185904456217153544375810280067230515294-x3_4, -326640557835176972818810921346435994618212511699701452909948165213793202867770210998489405471893729645622268115655639033758078636509576237302554335576631730085559843858617877502626166331583531556816922924864816637929054700592202533334089625644288273896168309864136047293358228290026459460082545647786035970941444350921016811485708567198207886588189896371078018397044806252246308227087658161951229143561271171716133962847995047566298735483509-x4_6, z_aux-1];
vars:=[x5_10, x5_9, x1_9, x5_8, x2_8, x1_8, x5_7, x4_7, x3_7, x2_7, x1_7, b_7, x5_6, x4_6, x3_6, x2_6, x1_6, b_6, x5_5, x4_5, x3_5, x2_5, x1_5, b_5, x5_4, x4_4, x3_4, x2_4, x1_4, b_4, x5_3, x4_3, x3_3, x2_3, x1_3, b_3, x5_2, x4_2, x3_2, x2_2, x1_2, b_2, x5_1, x4_1, x3_1, x2_1, x1_1, b_1, x5_0, x4_0, x3_0, x2_0, x1_0, b_0, z_aux, w_aux, a, d, k3, k4, k5, k6, k7];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [b, x3_4, x4_6] 100
# 171, [1 = 14, k5*x2_0 = 4, k6 = 4, x2_1*k5 = 3, k7*x1_5 = 2, k5*x2_2*x3_3 = 2, k5*x2_3*x3_2 = 2, x5_10 = 2, k5*x2_1*x4_5 = 2, k5*x2_2*x3_2 = 2, d*x2_3*k5 = 2, a*k5*x2_1 = 2, a*k5*x2_4 = 2, k6*x4_3 = 2, k5*x2_4*x3_2 = 2, x5_9 = 2, k6*x3_6 = 2, k4*x2_5 = 2, k5*x2_2*x3_1 = 2, k5*x2_0*x3_3 = 2, k5*x2_0*x4_1 = 2, k6*x4_0 = 2, k5*x2_0*x3_2 = 2, k5*x2_1*x3_5 = 2, x5_1 = 2, k5*x2_4*x4_0 = 2, k6*x3_3 = 2, k5*x2_2*x4_2 = 2, x1_5*k3 = 2, k4*x2_6 = 2, x1_6*k3 = 2, x5_4 = 2, k4*x2_7 = 2, k6*x4_1 = 2, k7*x1_6 = 2, k5*x2_5*x3_0 = 2, d*x2_4*k5 = 2, k5*x2_0*x3_5 = 2, k6*x4_4 = 2, k5*x2_0*x3_0 = 2, x1_1*k3 = 2, k5*x2_1*x4_3 = 2, k5*x2_6*x3_0 = 2, k7*x1_2 = 2, k4*x2_2 = 2, k5*x2_3*x3_0 = 2, x5_6 = 2, k5*x2_2*x4_1 = 2, k5*x2_6*x4_0 = 2, k5*x2_0*x4_3 = 2, k4*x2_3 = 2, k5*x2_4*x4_2 = 2, k5*x2_1*x3_3 = 2, k6*x3_1 = 2, k5*x2_1*x3_0 = 2, k5*x2_2*x4_4 = 2, k5*x2_0*x4_2 = 2, d*k5*x2_0 = 2, a*k5*x2_5 = 2, k5*x2_5*x4_0 = 2, k5*x2_3*x3_1 = 2, k5*x2_3*x4_2 = 2, k5*x2_1*x4_4 = 2, k6*x3_0 = 2, x5_8 = 2, k5*x2_3*x4_1 = 2, k5*x2_0*x4_0 = 2, k6*x3_2 = 2, x2_2*k5 = 2, k7*x1_7 = 2, k5*x2_0*x3_6 = 2, k4*x2_0 = 2, k5*x2_2*x3_0 = 2, k5*x2_4*x4_1 = 2, k7*x1_3 = 2, k5*x2_0*x4_5 = 2, k5*x2_1*x3_2 = 2, k5*x2_0*x3_1 = 2, k5*x2_2*x4_0 = 2, k5*x2_1*x4_0 = 2, k5*x2_1*x3_1 = 2, x1_4*k3 = 2, k4*x2_1 = 2, x2_5*d*k5 = 2, x1_1*k7 = 2, x2_6*d*k5 = 2, k3*x1_0 = 2, k6*x3_5 = 2, x5_2 = 2, x5_7 = 2, k5*x2_5*x3_1 = 2, k7*x1_0 = 2, k5*x2_0*x4_4 = 2, k5*x2_3*x3_3 = 2, k6*x4_2 = 2, k4*x2_4 = 2, k5*x2_4*x3_1 = 2, x5_5 = 2, x5_3 = 2, k6*x4_5 = 2, k5*x2_4*x3_0 = 2, k7*x1_8 = 2, a*k5*x2_6 = 2, k5*x2_2*x4_3 = 2, d*x2_1*k5 = 2, k5*x2_1*x4_2 = 2, x1_2*k3 = 2, a*k5*x2_3 = 2, a*k5*x2_2 = 2, k5*x2_5*x4_1 = 2, k5*x2_3*x4_0 = 2, k5*x2_1*x4_1 = 2, x1_7*k3 = 2, k5*x2_3*x4_3 = 2, k7*x1_4 = 2, a*k5*x2_0 = 2, d*x2_2*k5 = 2, x1_3*k3 = 2, x4_5 = 1, x2_3*k5 = 1, x2_2*k5*x4_5 = 1, x2_4 = 1, x4_1 = 1, x2_6*k5*x4_1 = 1, a*k5*x2_7 = 1, x2_6*k5*x3_1 = 1, k5*x2_0*x4_7 = 1, x1_3 = 1, x3_7 = 1, z_aux = 1, x5_0 = 1, x2_2 = 1, x2_8 = 1, x2_7*d*k5 = 1, x1_9 = 1, x2_5*k5*x3_2 = 1, x4_2 = 1, x2_2*k5*x3_5 = 1, x3_6 = 1, k6*x3_7 = 1, x4_4 = 1, x2_5*k5*x4_2 = 1, x1_2 = 1, x3_1 = 1, x2_6 = 1, k5*x2_7*x3_0 = 1, k7*x1_9 = 1, x4_3 = 1, k5*x2_7*x4_0 = 1, k4*x2_8 = 1, x4_7 = 1, x2_1 = 1, x2_4*k5*x3_3 = 1, x2_4*k5*x4_3 = 1, x3_3 = 1, x2_3 = 1, x3_5 = 1, x2_3*k5*x4_4 = 1, x2_5 = 1, x1_1 = 1, x1_8*k3 = 1, k6*x4_7 = 1, x2_1*k5*x3_6 = 1, x2_7 = 1, x1_4 = 1, k5*x2_0*x3_7 = 1, x1_7 = 1, x3_2 = 1, x1_6 = 1, x1_5 = 1, x1_8 = 1]
# 306, -5.039367613
# 118
# [b = 15, x4_6 = 6, x3_4 = 10]
# [b = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], x4_6 = [3, 2, 1, 3, 3, 2], x3_4 = [3, 2, 1, 3, 3, 2, 3, 3, 3, 3]]